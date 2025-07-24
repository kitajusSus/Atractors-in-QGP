# Krok 1: Importowanie potrzebnych pakietów
using DifferentialEquations
using Plots
using Random
using LaTeXStrings

# Używamy struktur, aby uniknąć zmiennych globalnych i elegancko grupować parametry.

"""
    HydroParams

Struktura przechowująca stałe fizyczne modelu hydrodynamiki (BRSSS).
"""
struct HydroParams
    C_τπ::Float64
    C_η::Float64
    C_λ1::Float64
end

"""
    SimSettings

Struktura przechowująca ustawienia dla przebiegu symulacji numerycznej.
"""
struct SimSettings
    n_points::Int
    τ_start::Float64
    τ_end::Float64
    tspan::Tuple{Float64, Float64}
    T_range::Tuple{Float64, Float64}
    A_range::Tuple{Float64, Float64}
end

# --- SEKCJA 2: FUNKCJE OBLICZENIOWE I WIZUALIZACYJNE ---

"""
    hydro_evolution!(du, u, p::HydroParams, τ)

Definicja układu równań różniczkowych dla modelu BRSSS.
Jest to "serce" modelu fizycznego.
"""
function hydro_evolution!(du, u, p::HydroParams, τ)
    T, A = u
    C_τπ, C_η, C_λ1 = p.C_τπ, p.C_η, p.C_λ1
    
    du[1] = (T / τ) * (-1/3 + A / 18)
    term_T = τ * T * (A + (C_λ1 / (12 * C_η)) * A^2)
    term_A2 = (2/9) * C_τπ * A^2
    du[2] = (1 / (C_τπ * τ)) * (8 * C_η - term_T - term_A2)
end

"""
    run_simulation(params::HydroParams, settings::SimSettings)

Uruchamia pełną symulację: generuje warunki początkowe, rozwiązuje
układ równań dla każdego z nich oraz dla trajektorii atraktora.

Zwraca: (solutions, sol_attractor)
"""
function run_simulation(params::HydroParams, settings::SimSettings)
    println("Rozpoczynanie symulacji z parametrami:")
    println("  C_τπ=$(round(params.C_τπ, digits=4)), C_η=$(round(params.C_η, digits=4)), C_λ1=$(round(params.C_λ1, digits=4))")
    println("  Przedział czasu τ: [$(settings.τ_start), $(settings.τ_end)]")

    Random.seed!(123)
    initial_conditions = [
        [rand(settings.T_range[1]:0.1:settings.T_range[2]), rand(settings.A_range[1]:0.1:settings.A_range[2])] 
        for _ in 1:settings.n_points
    ]

    solutions = []
    for u0 in initial_conditions
        prob = ODEProblem(hydro_evolution!, u0, settings.tspan, params)
        sol = solve(prob, Tsit5(), saveat=0.01)
        push!(solutions, sol)
    end

    # Wyznaczenie numerycznej trajektorii atraktora
    u0_attractor = [settings.T_range[2], 0.1] # Start blisko atraktora
    prob_attractor = ODEProblem(hydro_evolution!, u0_attractor, (0.01, settings.τ_end), params) 
    sol_attractor = solve(prob_attractor, Tsit5(), saveat=0.01)
    
    println("Zakończono rozwiązywanie równań.")
    return solutions, sol_attractor
end

"""
    create_log_ratio_animation(solutions, sol_attractor, settings, filename, fps)

Tworzy animację pokazującą logarytm stosunku anizotropii numerycznej do analitycznej.
"""
function create_log_ratio_animation(solutions, sol_attractor, settings, filename, fps)
    println("Generowanie animacji: $(filename)...")
    
    T_plot_min, T_plot_max = 50, 550
    τ_start, τ_end = settings.τ_start, settings.τ_end
    params = sol_attractor.prob.p # Pobranie parametrów z problemu
    C_η, C_τπ = params.C_η, params.C_τπ

    anim = @animate for τ_anim in τ_start:0.005:τ_end
        plot(
            xlabel="Temperatura T [MeV]",
            ylabel=L"ln(\frac{A_{num}}{A_{analit}})",
            title="Odległość od analitycznego atraktora w τ = $(round(τ_anim, digits=3)) fm/c",
            xlims=(T_plot_min, T_plot_max),
            ylims=(-0.5, 2.0),
            legend=:topright,
            grid=true, gridstyle=:dash, gridalpha=0.3
        )
        
        hline!([0], color=:magenta, linestyle=:dash, label="Atraktor Analityczny")
        
        for sol in solutions
            T_val, A_val = sol(τ_anim)
            A_attr = 8*C_η / (τ_anim * T_val) + (16*C_η * C_τπ) / (3 * (τ_anim * T_val)^2)
            scatter!([T_val], [log(A_val / A_attr)], 
                      marker_z=T_val, clims=(300, 500), 
                      markerstrokewidth=0, markersize=4, label="")
        end
    end

    gif(anim, filename, fps=fps)
    println("Animacja zapisana.")
end

"""
    create_log_distance_animation(solutions, sol_attractor, settings, filename, fps)

Tworzy animację pokazującą odległość od numerycznego atraktora w skali logarytmicznej.
"""
function create_log_distance_animation(solutions, sol_attractor, settings, filename, fps)
    println("Generowanie animacji: $(filename)...")

    τ_start, τ_end = settings.τ_start, settings.τ_end
    
    anim_log = @animate for τ_anim in τ_start:0.01:τ_end
        A_attr_num = sol_attractor(τ_anim)[2]
        
        plot(
            xlabel="Temperatura T [MeV]",
            ylabel="Odległość |A - A_attr|",
            title="Zbieżność do atraktora numerycznego w τ = $(round(τ_anim, digits=2)) fm/c",
            legend=false,
            yaxis=:log,
            ylims=(1e-4, 10.0)
        )
        
        hline!([1e-4], linestyle=:dash, color=:red)

        for sol in solutions
            T_val, A_val = sol(τ_anim)
            delta_A = abs(A_val - A_attr_num) + 1e-9
            initial_A = sol.u[1][2]
            scatter!([T_val], [delta_A], marker_z=initial_A, clims=settings.A_range, markerstrokewidth=0, markersize=4)
        end
    end

    gif(anim_log, filename, fps=fps)
    println("Animacja zapisana.")
end


# --- SEKCJA 3: GŁÓWNA FUNKCJA STERUJĄCA ---

"""
    main()

Główna funkcja programu. Definiuje parametry, uruchamia symulację i wizualizacje.
Idealne miejsce do eksperymentowania.
"""
function main()
    println("--- Symulacja Hydrodynamicznego Atraktora ---")

    # === MIEJSCE DO EKSPERYMENTOWANIA ===
    # Zmień parametry tutaj, aby zobaczyć, jak wpływają na wyniki.
    
    # Zestaw parametrów fizycznych (zgodnie z Teorią N=4 SYM)
    hydro_params = HydroParams(
        (2 - log(2)) / (2 * π),  # C_τπ
        1 / (4 * π),            # C_η
        1 / (2 * π)             # C_λ1
    )

    # Ustawienia symulacji
    sim_settings = SimSettings(
        200,                    # n_points
        0.2,                    # τ_start
        1.0,                    # τ_end
        (0.2, 1.0),             # tspan
        (300.0, 500.0),         # T_range
        (0.0, 4.0)              # A_range
    )
    # ====================================

    # Uruchomienie obliczeń
    solutions, sol_attractor = run_simulation(hydro_params, sim_settings)

    # Tworzenie wizualizacji
    create_log_ratio_animation(
        solutions, sol_attractor, sim_settings, 
        "anim_log_ratio.gif", 15
    )
    
    create_log_distance_animation(
        solutions, sol_attractor, sim_settings, 
        "anim_log_distance.gif", 15
    )

    println("--- Zakończono pomyślnie ---")
end
