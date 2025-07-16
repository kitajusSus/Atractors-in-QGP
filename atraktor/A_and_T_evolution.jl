# 1. Importowanie potrzebnych pakietów
using DifferentialEquations
using Plots
using Printf

# 2. Definicja stałych fizycznych
# Używamy wartości dla teorii N=4 SYM jako przykładu (z artykułu przeglądowego, równanie 33)
# Zapewniają one realistyczne zachowanie atraktora.
const C_η = 1 / (4π)
const C_τπ = (2 - log(2)) / (2π)
const C_λ1 = 1 / (2π)
const PARAMS = (C_η, C_τπ, C_λ1)

# 3. Definicja układu równań różniczkowych (ODE)
# u = [T, A] - wektor stanu (Temperatura, Anizotropia)
# p - parametry fizyczne
# τ - czas właściwy
function hydro_evolution!(du, u, p, τ)
    T, A = u
    C_η, C_τπ, C_λ1 = p

    # Równanie na d(T)/dτ (z τ d(logT)/dτ = -1/3 + A/18)
    du[1] = (T / τ) * (-1/3 + A / 18)

    # Równanie na d(A)/dτ (z teorii MIS)
    # Zmieniona postać równania (34) z artykułu
    term_T_dep = τ * T * (A + (C_λ1 / (12 * C_η)) * A^2)
    term_A2_dep = (2/9) * C_τπ * A^2
    du[2] = (1 / (C_τπ * τ)) * (8*C_η - term_T_dep - term_A2_dep)
end

# 4. Parametry symulacji
const N_POINTS = 500         # Liczba generowanych warunków początkowych
const T_RANGE = (200.0, 500.0) # Zakres początkowej temperatury T [MeV]
const A_RANGE = (-1.0, 5.0)   # Zakres początkowej anizotropii A

const τ_start = 0.2           # Czas początkowy [fm/c]
const τ_end = 1.2             # Czas końcowy [fm/c]
const τ_span = (τ_start, τ_end)

# Momenty czasu, dla których będziemy tworzyć wykresy ("migawki")
const τ_snapshots = [0.2, 0.4, 0.6, 0.8, 1.0, 1.2]

# 5. Generowanie warunków początkowych i rozwiązywanie ODE
println("Rozpoczynam symulację ewolucji dla $N_POINTS warunków początkowych...")

# Słownik do przechowywania wyników dla każdej migawki czasu
# Klucz: czas τ, Wartość: (wektor T, wektor A, wektor T_początkowych do kolorowania)
results = Dict{Float64, Tuple{Vector{Float64}, Vector{Float64}, Vector{Float64}}}()
for τ_val in τ_snapshots
    results[τ_val] = ([], [], [])
end

# Generujemy losowe warunki początkowe
initial_conditions_T = T_RANGE[1] .+ (T_RANGE[2] - T_RANGE[1]) .* rand(N_POINTS)
initial_conditions_A = A_RANGE[1] .+ (A_RANGE[2] - A_RANGE[1]) .* rand(N_POINTS)

# Pętla po wszystkich warunkach początkowych
for i in 1:N_POINTS
    T₀ = initial_conditions_T[i]
    A₀ = initial_conditions_A[i]
    u₀ = [T₀, A₀]

    # Definiowanie i rozwiązywanie problemu ODE
    prob = ODEProblem(hydro_evolution!, u₀, τ_span, PARAMS)
    
    # Rozwiązujemy równanie, mówiąc solverowi, żeby zapisał wyniki w punktach τ_snapshots
    sol = solve(prob, Tsit5(), saveat=τ_snapshots)

    # Zapisujemy wyniki dla każdej migawki
    for (j, τ_val) in enumerate(τ_snapshots)
        T_at_τ, A_at_τ = sol.u[j]
        push!(results[τ_val][1], T_at_τ)
        push!(results[τ_val][2], A_at_τ)
        push!(results[τ_val][3], T₀) # Zapisujemy T₀ do kolorowania
    end
    
    # Wyświetlanie postępu
    if i % 50 == 0
        println("... Obliczono $i / $N_POINTS trajektorii.")
    end
end

println("Symulacja zakończona. Generowanie animacji...")

# 6. Tworzenie animacji z serii wykresów
anim = @animate for τ_val in τ_snapshots
    
    # Pobieranie danych dla danego czasu
    T_vals, A_vals, T0_vals = results[τ_val]

    # Tworzenie pojedynczego wykresu (jednej klatki animacji)
    scatter(T_vals, A_vals,
        title = @sprintf("Przestrzeń Stanów (T, A) w τ = %.2f fm/c", τ_val),
        xlabel = "Temperatura T [MeV]",
        ylabel = "Anizotropia A",
        label = "",
        xlims = (0, T_RANGE[2]), # Stałe granice osi dla lepszego porównania
        ylims = (A_RANGE[1], A_RANGE[2]),
        markersize = 3,
        markerstrokewidth = 0,
        alpha = 0.7,
        zcolor = T0_vals, # Kolorowanie punktów wg ich początkowej temperatury
        colorbar_title = " T₀ [MeV]",
        framestyle = :box
    )
end

# Zapisywanie animacji do pliku GIF
gif(anim, "ewolucja_atraktora.gif", fps = 1)

println("\nAnimacja została zapisana do pliku 'ewolucja_atraktora.gif'")

# Opcjonalnie: można też wygenerować osobne pliki dla każdego wykresu
# for τ_val in τ_snapshots
#     T_vals, A_vals, T0_vals = results[τ_val]
#     p = scatter(T_vals, A_vals, ..., ) # jak wyżej
#     savefig(p, "stan_w_tau_$(replace(string(τ_val), "." => "_")).png")
# end
