# Plik: HydroSim.jl
module modHydroSim

# --- Zależności ---
using DifferentialEquations
using Plots
using Random
using LaTeXStrings

# --- Publiczny Interfejs Modułu ---
export HydroParams, SimSettings, SimResult,
       PARAMS_SYM_THEORY, PARAMS_MIS_TOY_MODEL,
       run_simulation, 
       create_log_ratio_animation, create_log_distance_animation


"""
    HydroParams

Niezmienna struktura przechowująca stałe fizyczne modelu hydrodynamiki (BRSSS).
Jej pól nie można zmienić po utworzeniu obiektu.
"""
struct HydroParams
    C_τπ::Float64
    C_η::Float64
    C_λ1::Float64
end

# --- Predefiniowane, niezmienne zestawy parametrów ---
"Parametry z teorii N=4 SYM (holografia/AdS-CFT), pełny model BRSSS."
const PARAMS_SYM_THEORY = HydroParams(
    (2 - log(2)) / (2 * π),  # C_τπ
    1 / (4 * π),            # C_η
    1 / (2 * π)             # C_λ1
)

"Parametry z pracy PRL (Heller, Spaliński et al.), model zabawkowy MIS."
const PARAMS_TOY_MODEL = HydroParams(1.0, 0.75, 0.0)
###############################################################
"""
    SimSettings

Niezmienna struktura przechowująca ustawienia symulacji.
Posiada konstruktor z argumentami kluczowymi dla łatwego tworzenia.
"""
struct SimSettings
    n_points::Int
    τ_start::Float64
    τ_end::Float64
    tspan::Tuple{Float64, Float64}
    T_range::Tuple{Float64, Float64}
    A_range::Tuple{Float64, Float64}
end
# Konstruktor z argumentami kluczowymi i wartościami domyślnymi
function SimSettings( 
        n_points=200, τ_start=0.2, τ_end=1.2,
        T_range=(300.0, 500.0), A_range=(0.0, 4.0))
    return settings = (n_points, τ_start, τ_end, (τ_start, τ_end), T_range, A_range)
end


"""
    SimResult

Struktura przechowująca kompletne wyniki symulacji wraz z metadanymi.
"""
struct SimResult
    solutions::Vector{ODESolution}
    sol_attractor::ODESolution
    params::HydroParams
    settings::SimSettings
end

# Funkcja wewnętrzna, nieeksportowana
function _hydro_evolution!(du, u, p::HydroParams, τ)
    T, A = u
    C_τπ, C_η, C_λ1 = p.C_τπ, p.C_η, p.C_λ1
    du[1] = (T / τ) * (-1/3 + A / 18)
    term_T = τ * T * (A + (C_λ1 / (12 * C_η)) * A^2)
    term_A2 = (2/9) * C_τπ * A^2
    du[2] = (1 / (C_τπ * τ)) * (8 * C_η - term_T - term_A2)
end

"""
    run_simulation(params::HydroParams, settings::SimSettings) -> SimResult

Uruchamia pełną symulację i zwraca jeden obiekt z wynikami.
"""
function run_simulation(params::HydroParams, settings::SimSettings)
    println("--- Rozpoczynanie Obliczeń Numerycznych...")
    rng = Xoshiro(5)
    initial_conditions = [
        [rand(rng, settings.T_range[1]:0.1:settings.T_range[2]), rand(rng, settings.A_range[1]:0.1:settings.A_range[2])] 
        for _ in 1:settings.n_points
    ]

    solutions = ODESolution[]
    for u0 in initial_conditions
        prob = ODEProblem(_hydro_evolution!, u0, settings.tspan, params)
        sol = solve(prob, Tsit5(), saveat=0.01)
        push!(solutions, sol)
    end

    u0_attractor = [settings.T_range[2], 0.1]
    prob_attractor = ODEProblem(_hydro_evolution!, u0_attractor, (0.01, settings.τ_end), params) 
    sol_attractor = solve(prob_attractor, Tsit5(), saveat=0.01)
    
    println("--- Obliczenia Zakończone. ---")
    return SimResult(solutions, sol_attractor, params, settings)
end

"""
    create_log_ratio_animation(result::SimResult; ...)

Tworzy animację stosunku A_num / A_analit.
"""
function create_log_ratio_animation(result::SimResult; filename="anim_log_ratio.gif", fps=15)
    println("Generowanie animacji: $(filename)...")
    
    C_η, C_τπ = result.params.C_η, result.params.C_τπ

    anim = @animate for τ_anim in result.settings.τ_start:0.005:result.settings.τ_end
        plot(
            xlabel="Temperatura T [MeV]", ylabel=L"ln(A_{num} / A_{analit})",
            title="Odległość od analitycznego atraktora w τ = $(round(τ_anim, digits=3)) fm/c",
            xlims=(50, 850), ylims=(-0.5, 2.0), legend=:topright,
            grid=true, gridstyle=:dash, gridalpha=0.3
        )
        hline!([0], color=:magenta, linestyle=:dash, label="Atraktor Analityczny")
        
        for sol in result.solutions
            T_val, A_val = sol(τ_anim)
            A_attr_analit = 8*C_η / (τ_anim * T_val) + (16*C_η * C_τπ) / (3 * (τ_anim * T_val)^2)
            scatter!([T_val], [log(A_val / A_attr_analit)], 
                      marker_z=T_val, clims=result.settings.T_range, 
                      markerstrokewidth=0, markersize=4, label="")
        end
    end

    gif(anim, filename, fps=fps)
    println("Animacja zapisana: $(filename)")
end
"""
    create_log_ratio_animation(solutions, sol_attractor, settings; filename="anim_log_ratio.gif", fps=15)

Tworzy animację stosunku `A_num/A_analit`. Funkcja jest **szybka**,
ponieważ operuje na gotowych danych.
"""
function create_log_ratio_animation(solutions, settings; filename="anim_log_ratio.gif", fps=15)
    println("Generowanie animacji: $(filename)...")
    
    params = solutions[1].prob.p
    C_η, C_τπ = params.C_η, params.C_τπ

    anim = @animate for τ_anim in settings.τ_start:0.005:settings.τ_end
        plot(
            xlabel="Temperatura T [MeV]",
            ylabel=L"ln(A_{num} / A_{analit})",
            title="Odległość od analitycznego atraktora w τ = $(round(τ_anim, digits=3)) fm/c",
            xlims=(50, 550), ylims=(-0.5, 2.0), legend=:topright,
            grid=true, gridstyle=:dash, gridalpha=0.3
        )
        hline!([0], color=:magenta, linestyle=:dash, label="Atraktor Analityczny")
        
        for sol in solutions
            T_val, A_val = sol(τ_anim)
            A_attr = 8*C_η / (τ_anim * T_val) + (16*C_η * C_τπ) / (3 * (τ_anim * T_val)^2)
            scatter!([T_val], [log(A_val / A_attr)], 
                      marker_z=T_val, clims=settings.T_range, 
                      markerstrokewidth=0, markersize=4, label="")
        end
    end

    gif(anim, filename, fps=fps)
    println("Animacja zapisana: $(filename)")
end

"""
    create_log_distance_animation(solutions, sol_attractor, settings; filename="anim_log_distance.gif", fps=15)

Tworzy animację odległości od numerycznego atraktora. Funkcja jest **szybka**.
"""
function create_log_distance_animation(solutions, sol_attractor, settings; filename="anim_log_distance.gif", fps=15)
    println("Generowanie animacji: $(filename)...")
    
    anim = @animate for τ_anim in settings.τ_start:0.01:settings.τ_end
        A_attr_num = sol_attractor(τ_anim)[2]
        
        plot(
            xlabel="Temperatura T [MeV]", ylabel="Odległość |A - A_attr|",
            title="Zbieżność do atraktora numerycznego w τ = $(round(τ_anim, digits=2)) fm/c",
            legend=false, yaxis=:log, ylims=(1e-4, 10.0)
        )
        hline!([1e-4], linestyle=:dash, color=:red)

        for sol in solutions
            T_val, A_val = sol(τ_anim)
            delta_A = abs(A_val - A_attr_num) + 1e-9
            scatter!([T_val], [delta_A], marker_z=sol.u[1][2], clims=settings.A_range, 
                      markerstrokewidth=0, markersize=4)
        end
    end

    gif(anim, filename, fps=fps)
    println("Animacja zapisana: $(filename)")
end

end
