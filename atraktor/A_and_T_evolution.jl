# --- Załadowanie potrzebnych pakietów ---
using DifferentialEquations
using Plots
using Printf

# Ustawienie domyślnego backendu na GR (szybki)
gr()

# --- Definicje fizyczne ---
const Cη = 1 / (4π)
const Cτ = (2 - log(2)) / (2π)
const ħc = 197.3  # Stała konwersji MeV·fm

# --- Układ równań różniczkowych (ewolucja w czasie τ) ---
# To jest serce naszej symulacji.
# u[1] -> A (anizotropia), u[2] -> T (temperatura w MeV)
function attractor_system!(du, u, p, τ)
    A, T = u

    # Zmienna pomocnicza w. Zabezpieczenie przed dzieleniem przez zero.
    w = τ * T / ħc
    if w < 1e-12
        w = 1e-12
    end

    # Równanie na dT/dτ
    # Zabezpieczenie przed A ≈ -12
    denom_T = τ * (A + 12.0)
    if abs(denom_T) < 1e-12
        denom_T = sign(denom_T) * 1e-12
    end
    dT_dτ = T * (A - 6.0) / denom_T
    du[2] = dT_dτ

    # Równanie na dA/dτ (z reguły łańcuchowej)
    # 1. Obliczamy dA/dw
    denom_A_dw = Cτ * (1.0 + A / 12.0)
    if abs(denom_A_dw) < 1e-12
        denom_A_dw = sign(denom_A_dw) * 1e-12
    end
    num_A_dw = (3.0 / 2.0) * (8.0 * Cη / w - A) - (Cτ / (3.0 * w)) * A^2
    dA_dw = num_A_dw / denom_A_dw

    # 2. Obliczamy dw/dτ
    dw_dτ = (T + τ * dT_dτ) / ħc
    
    # 3. Łączymy wszystko
    du[1] = dA_dw * dw_dτ
end

# --- Parametry symulacji ---
const N_POINTS = 10     # Liczba symulowanych trajektorii
const τ_initial = 0.2      # Początkowy czas fizyczny [fm/c]
const τ_final = 1.5        # Końcowy czas fizyczny [fm/c]

# Zakresy dla losowych warunków początkowych
const T_range = (150.0, 400.0) # [MeV]
const A_range = (-2.0, 5.0)

# --- Pętla rozwiązująca ---
println("Generowanie i rozwiązywanie $N_POINTS trajektorii...")
solutions = []
for i in 1:N_POINTS
    print("Rozwiązywanie trajektorii $i/$N_POINTS\r")

    # Losujemy stan początkowy (T₀, A₀) w czasie τ_initial
    T₀ = rand(T_range)
    A₀ = rand(A_range)
    u₀ = [A₀, T₀]

    # Definiujemy problem do rozwiązania
    τ_span = (τ_initial, τ_final)
    prob = ODEProblem(attractor_system!, u₀, τ_span)

    # Rozwiązujemy, używając odpornego solvera
    # `save_everystep=false` i `dense=true` są kluczowe dla wydajności
    # i dokładnej interpolacji później.
    sol = solve(prob, Rodas5(autodiff=false), save_everystep=false, dense=true, abstol=1e-6, reltol=1e-6, maxiters=1e7)

    # Dodajemy tylko pomyślnie zakończone rozwiązania
    if sol.retcode == :Success
        push!(solutions, sol)
    end
end
println("\nZakończono. Pomyślnie rozwiązano $(length(solutions)) trajektorii.")

# --- Tworzenie animacji ---
if isempty(solutions)
    error("Nie udało się rozwiązać żadnej trajektorii. Nie można stworzyć animacji.")
end

# Definiujemy klatki czasowe dla naszej animacji
τ_frames = range(τ_initial, τ_final, length=150)

println("Tworzenie animacji z $(length(τ_frames)) klatek...")
anim = @animate for (i, τ_target) in enumerate(τ_frames)
    print("Generowanie klatki $i/$(length(τ_frames))\r")

    # Zbieramy współrzędne wszystkich punktów w tej klatce
    points_T = Float64[]
    points_A = Float64[]

    for sol in solutions
        # Sprawdzamy, czy ten moment czasowy jest w zakresie danego rozwiązania
        if τ_initial <= τ_target <= sol.t[end]
            # Używamy interpolacji `sol(τ_target)`, aby uzyskać dokładne wartości
            A_target, T_target = sol(τ_target)
            
            # Dodajemy tylko sensowne punkty do narysowania
            if isfinite(A_target) && isfinite(T_target)
                push!(points_A, A_target)
                push!(points_T, T_target)
            end
        end
    end

    # Rysujemy klatkę (za każdym razem od nowa)
    plot(
        size=(800, 600),
        title="Przestrzeń Stanów (T, A) w τ = $(@sprintf("%.2f", τ_target)) fm/c",
        xlabel="Temperatura T [MeV]",
        ylabel="Anizotropia A",
        xlims=(0, 500),
        ylims=(-2, 5),
        legend=false,
        grid=true
    )

    # Rysujemy punkty, jeśli jakieś mamy
    if !isempty(points_T)
        scatter!(
            points_T, 
            points_A, 
            marker_z=points_T, # Koloruj punkty wg temperatury
            color=:viridis,
            clims=T_range,
            markersize=4,
            markerstrokewidth=0
        )
    end
end

# Zapisanie gotowej animacji jako plik GIF
gif_filename = "attractor_evolution_final_clean.gif"
gif(anim, gif_filename, fps=20)

println("\nAnimacja została pomyślnie zapisana jako: $gif_filename")
