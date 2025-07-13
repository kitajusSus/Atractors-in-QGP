# -----------------------------------------------------------------
# SEKCJA 1: GENEROWANIE DANYCH - SYMULACJA HYDRODYNAMICZNA
# -----------------------------------------------------------------
using DifferentialEquations, Plots, MultivariateStats, Statistics, Flux

# Definicja równania i parametrów
const C_eta = 1 / (4π)
const C_tau_pi = (2 - log(2)) / (2π)
const C_lambda_1 = 1 / (2π)
const p_hydro = (C_eta, C_tau_pi, C_lambda_1)

# Uproszczone równanie różniczkowe (model MIS, C_lambda_1 = 0)
# jest bardziej stabilne numerycznie i wystarczające do demonstracji.
function mis_flow_simple!(du, u, p, t)
    T, T_dot = u
    C_eta, C_tau_pi, _ = p # Ignorujemy C_lambda_1

    T_ddot = -( (3/2)/t * (2*C_tau_pi) * T_dot^2 +
                ( (T/C_eta) * C_eta + 11*C_tau_pi/(3*t) ) * T_dot +
                ( T^2 / (3*C_eta*t) - 4*(C_eta-C_tau_pi)*T/(9*C_tau_pi*t) ) )
    
    du[1] = T_dot
    du[2] = T_ddot
end

# Generowanie "chmury" punktów - trajektorii
τ_span = (0.5, 5.0)
num_trajectories = 20000
T0 = 1 # GeV

# !!! KLUCZOWA ZMIANA: Znacznie zawężony i bardziej fizyczny zakres warunków początkowych !!!
T_dot0_range = -10:0.005:10
initial_conditions = [[T0, T_dot0] for T_dot0 in rand(T_dot0_range, num_trajectories)]

# Zbieramy wszystkie punkty w jednej macierzy dla każdej "migawki" czasowej
τ_snapshots = 0.5:0.5:5.0
all_snapshots = []

println("Generowanie danych z symulacji hydrodynamicznej...")
successful_trajectories = 0
for ic in initial_conditions
    prob = ODEProblem(mis_flow_simple!, ic, τ_span, p_hydro)
    # Zwiększamy tolerancję, żeby więcej solverów się powiodło
    sol = solve(prob, Tsit5(), saveat=τ_snapshots, abstol=1e-6, reltol=1e-6)
    if sol.retcode == :Success
        push!(all_snapshots, hcat(sol.u...))
        global successful_trajectories += 1
    end
end
println("Dane wygenerowane. Liczba udanych trajektorii: $successful_trajectories")

# !!! OBSŁUGA BŁĘDU: Sprawdzamy, czy mamy jakiekolwiek dane !!!
if successful_trajectories == 0
    error("Nie udało się wygenerować żadnej stabilnej trajektorii. Spróbuj jeszcze bardziej zawęzić T_dot0_range lub zwiększyć τ0.")
end

data_snapshots = [hcat([snapshot[:,i] for snapshot in all_snapshots]...) for i in 1:length(τ_snapshots)]

# -----------------------------------------------------------------
# SEKCJA 2: ANALIZA PCA
# -----------------------------------------------------------------
println("\nPrzeprowadzanie analizy PCA...")
pca_explained_variances = []

for snapshot in data_snapshots
    # !!! POPRAWKA: `mean=0` zostanie automatycznie zinterpretowane poprawnie,
    # ale można jawnie podać wektor zer.
    M = fit(PCA, snapshot; pratio=1.0, mean=mean(snapshot, dims=2))
    push!(pca_explained_variances, principalvars(M) ./ tvar(M))
end

# Wykres ewolucji wariancji wyjaśnionej
plot_pca = plot(τ_snapshots, [ev[1] for ev in pca_explained_variances],
    label="PCA: Składowa główna 1", lw=3,
    xlabel="Czas własny τ [fm/c]", ylabel="Wariancja wyjaśniona",
    title="Dynamiczna Redukcja Wymiarowości (PCA vs Autoenkoder)",
    legend=:right, framestyle=:box, ylims=(0, 1.05)
)
plot!(plot_pca, τ_snapshots, [ev[2] for ev in pca_explained_variances],
    label="PCA: Składowa główna 2", lw=3
)
println("Analiza PCA zakończona.")

# -----------------------------------------------------------------
# SEKCJA 3: BUDOWA I TRENING AUTOENKODERA
# -----------------------------------------------------------------
println("\nBudowa i trening autoenkodera...")

# Definicja architektury Autoenkodera
latent_dim = 1
# Normalizacja danych przed podaniem do sieci neuronowej jest kluczowa!
struct Normalizer
    μ::Matrix{Float32}
    σ::Matrix{Float32}
end
(n::Normalizer)(x) = (x .- n.μ) ./ (n.σ .+ 1f-8) # Dodajemy epsilon, by uniknąć dzielenia przez 0
(n::Normalizer)(x, rev::Bool) = rev ? (x .* n.σ) .+ n.μ : n(x)

# Użyjemy mniejszej sieci, co powinno wystarczyć i przyspieszyć trening
encoder = Chain(Dense(2 => 8, relu), Dense(8 => latent_dim))
decoder = Chain(Dense(latent_dim => 8, relu), Dense(8 => 2))

# Funkcja straty
loss(x, norm) = Flux.mse(decoder(encoder(norm(x))), norm(x))

# Pętla treningowa
ae_reconstruction_errors = []
for snapshot in data_snapshots
    # Normalizacja danych dla każdej migawki
    μ = mean(snapshot, dims=2)
    σ = std(snapshot, dims=2)
    norm = Normalizer(μ, σ)
    
    # Resetowanie parametrów modelu dla każdej nowej migawki
    global encoder = Chain(Dense(2 => 8, relu), Dense(8 => latent_dim))
    global decoder = Chain(Dense(latent_dim => 8, relu), Dense(8 => 2))
    
    opt = ADAM(0.001)
    train_data = [(Float32.(snapshot),)] # Dane muszą być odpowiedniego typu

    # Trenujemy autoenkoder
    Flux.@epochs 200 Flux.train!((x) -> loss(x, norm), Flux.params(encoder, decoder), train_data, opt)

    reconstruction_error = loss(Float32.(snapshot), norm)
    # Obliczamy całkowitą wariancję na znormalizowanych danych
    total_variance = sum(var(norm(snapshot), dims=2))
    
    push!(ae_reconstruction_errors, reconstruction_error / total_variance)
end

ae_explained_variances = 1.0 .- ae_reconstruction_errors

plot!(plot_pca, τ_snapshots, ae_explained_variances,
    label="Autoenkoder (1D latent)", lw=3, color=:red, linestyle=:dash
)
println("Trening autoenkodera zakończony.")

# -----------------------------------------------------------------
# SEKCJA 4: WIZUALIZACJA I PORÓWNANIE
# -----------------------------------------------------------------

display(plot_pca)

println("\nGenerowanie wizualizacji zapadania się chmury punktów...")
anim = @animate for i in 1:length(τ_snapshots)
    scatter(data_snapshots[i][1,:], data_snapshots[i][2,:],
        xlabel="Temperatura T [GeV]",
        ylabel="Pochodna Ṫ [GeV/fm]",
        title="Ewolucja w przestrzeni fazowej (τ = $(round(τ_snapshots[i], digits=1)) fm/c)",
        xlims=(0.1, 0.7), ylims=(-0.5, 0.5), # Dostosowane do nowych warunków początkowych
        legend=false, framestyle=:box,
        marker_z=1:successful_trajectories, color=:viridis, markerstrokewidth=0
    )
end

gif(anim, "atraktor_animacja.gif", fps=2)
println("Animacja 'atraktor_animacja.gif' została zapisana.")
