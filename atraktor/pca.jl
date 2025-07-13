# =================================================================
# PROJEKT: ODKRYWANIE NIELINIOWEJ GEOMETRII ATRAKTORA (Wersja 2.2 - Poprawne API Flux)
# Język: Julia
# Pakiety: Plots, Flux, MultivariateStats, Random, Statistics
# =================================================================

using Plots, Flux, MultivariateStats, Random, Statistics

# Ustawienie globalnego stylu wykresów dla spójności
theme(:default)
const C_eta = 1 / (4π)
const C_tau_pi = (2 - log(2)) / (2π)
const C_lambda_1 = 1 / (2π)

w_true = range(0.1, 10, length=1000)
# brsss
A_true = 8 .* C_eta ./ w_true .+ (16 .* C_eta .* (C_tau_pi - C_lambda_1)) ./ (3 .* w_true.^2)
# -----------------------------------------------------------------
# KROK 1: GENERACJA DANYCH TESTOWYCH
# -----------------------------------------------------------------
println("Generowanie danych testowych...")

function generate_data(num_points=1000, noise_level=0.5)




    #A_true = 1.2 ./ w_true .+ 0.3 * sin.(2 .* w_true)
    w_data, A_data = Float32[], Float32[]

    for _ in 1:num_points
        idx = rand(1:length(w_true))
        w_base, A_base = w_true[idx], A_true[idx]
        current_noise = randn() * noise_level / (1 + w_base / 2)
        #dodawanie szumu do tablic danych float32
        push!(w_data, w_base + current_noise * 0.2)
        push!(A_data, A_base + current_noise)
    end
    
    return hcat(w_data, A_data)'
end

data_cloud = generate_data()
println("Dane wygenerowane. Wymiary: $(size(data_cloud))")

# Wizualizacja 1: Wygenerowane dane
p1 = scatter(data_cloud[1,:], data_cloud[2,:],
    label="Dane 'z symulacji'", markersize=2, markerstrokewidth=0, alpha=0.6,
    xlabel="Bezwymiarowy czas w", ylabel="Anizotropia A",
    title="1. Wygenerowana chmura danych", legend=:topright
)
display(p1)
#println(" Generowanie danych zakonczone. (zobacz wykres) Naciśnij Enter, aby kontynuować...")
#readline()


# -----------------------------------------------------------------
# KROK 2: ANALIZA PCA (PODEJŚCIE LINIOWE)
# -----------------------------------------------------------------
println("\nPrzeprowadzanie analizy PCA...")
M = fit(PCA, data_cloud; maxoutdim=1)
data_reconstructed_pca = reconstruct(M, transform(M, data_cloud))
mse_pca = mean((data_cloud .- data_reconstructed_pca).^2)
println("Analiza PCA zakończona. Błąd rekonstrukcji: $mse_pca")

p2 = scatter(data_cloud[1,:], data_cloud[2,:],
    label="Dane oryginalne", color=:gray, alpha=0.3, markersize=2)
scatter!(p2, data_reconstructed_pca[1,:], data_reconstructed_pca[2,:],
    label="Rekonstrukcja PCA (Liniowy atraktor)", color=:blue, markersize=2,
    xlabel="w", ylabel="A", title="2. Liniowy Atraktror znaleziony przez PCA"
)
display(p2)
#println("analiza PCA skonczona. (wykres) Naciśnij Enter, aby kontynuować...")
#readline()

# -----------------------------------------------------------------
# KROK 3: AUTOENKODER (PODEJŚCIE NIELINIOWE)
# -----------------------------------------------------------------
println("\nBudowa i trening autoenkodera ( biblioteka Flux)...")

# Normalizacja danych
μ = mean(data_cloud, dims=2)
σ = std(data_cloud, dims=2)
data_normalized = Float32.((data_cloud .- μ) ./ σ)

# Architektura sieci
latent_dim = 1
autoencoder = Chain(
    Dense(2 => 16, tanh),
    Dense(16 => 8, tanh),
    Dense(8 => latent_dim),
    Dense(latent_dim => 8, tanh),
    Dense(8 => 16, tanh),
    Dense(16 => 2)
)

# Definicja funkcji straty 
# Teraz przyjmuje 3 argumenty: model, wejście (x) i oczekiwane wyjście (y)
# W autoenkoderze x i y to te same dane.
loss(model, x, y) = Flux.mse(model(x), y)

# Optymalizator
opt_state = Flux.setup(Adam(0.001), autoencoder)

# Przygotowanie danych do treningu
# DataLoader jest lepszy do zarządzania danymi, zwłaszcza większymi
train_loader = Flux.DataLoader((data_normalized, data_normalized), batchsize=128, shuffle=true)

# Pętla treningowa
println("Trening...")
epochs = 100
anim_training = @animate for i in 1:epochs
    #println("Epoka $i / $epochs")
    # Pętla po minibatchach danych
    for (x_batch, y_batch) in train_loader
        Flux.train!(loss, autoencoder, [(x_batch, y_batch)], opt_state)
    end
    
    if (i % 20 == 0)
        reconstructed_norm = autoencoder(data_normalized)
        reconstructed_denorm = (reconstructed_norm .* σ) .+ μ
        
        p_train = scatter(data_cloud[1,:], data_cloud[2,:], color=:gray, alpha=0.2, label="", framestyle=:box)
        scatter!(p_train, reconstructed_denorm[1,:], reconstructed_denorm[2,:],
            label="Rekonstrukcja AE", color=:magenta, markersize=2.5,
            xlabel="w", ylabel="A", title="3. Trening Autoenkodera (Epoka $i / $epochs)",
            xlims=extrema(data_cloud[1,:]), ylims=extrema(data_cloud[2,:])
        )
    end
end

gif(anim_training, "trening_autoenkodera.gif", fps=15)
println("Animacja treningu zapisana jako 'trening_autoenkodera.gif'.")

# Finalna rekonstrukcja i błąd
data_reconstructed_ae = (autoencoder(data_normalized) .* σ) .+ μ
mse_ae = mean((data_cloud .- data_reconstructed_ae).^2)
println("Trening zakończony. Błąd rekonstrukcji: $mse_ae")

# -----------------------------------------------------------------
# KROK 4: WIZUALIZACJA KOŃCOWA I PORÓWNANIE
# -----------------------------------------------------------------
p_final = scatter(data_cloud[1,:], data_cloud[2,:],
    label="Dane 'z symulacji'",
    markersize=2, markerstrokewidth=0, alpha=0.3, color=:gray,
    xlabel="Bezwymiarowy czas w", ylabel="Anizotropia A",
    title="Porównanie: PCA vs Autoenkoder",
    legend=:topright, framestyle=:box
)
scatter!(p_final, data_reconstructed_pca[1,:], data_reconstructed_pca[2,:],
    label="PCA (MSE: $(round(mse_pca, digits=4)))",
    markersize=2.5, color=:blue
)
scatter!(p_final, data_reconstructed_ae[1,:], data_reconstructed_ae[2,:],
    label="Autoenkoder (MSE: $(round(mse_ae, digits=4)))",
    markersize=2.5, color=:magenta
)
plot!(p_final, range(0.1, 10, length=1000), A_true, 
label="A_true (bez szumu)", color=:red, linewidth=2)

display(p_final)
savefig(p_final, "wykresy_pca_ae.png")

println("\n================ WNIOSKI ================")
if mse_ae < mse_pca * 0.5
    println(">>> Autoenkoder osiągnął co najmniej 2x niższy błąd rekonstrukcji.")
    println(">>> Sugeruje to, że ukryty atraktor w danych jest strukturą NIELINIOWĄ.")
    println(">>> Liniowe metody, jak PCA, nie są w stanie w pełni uchwycić jego geometrii.")
else
    println(">>> Błędy obu metod są porównywalne.")
    println(">>> Wskazuje to, że atraktor w tym zbiorze danych jest w przybliżeniu liniowy,")
    println(">>> lub zastosowana sieć neuronowa jest zbyt prosta, aby nauczyć się nieliniowości.")
end
println("=========================================")
