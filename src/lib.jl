# Plik: lib.jl
module modHydroSim

# --- Zależności ---
using DifferentialEquations
using Random
using Distributions
using Plots


# --- Publiczny Interfejs Modułu ---
export HydroParams, SimSettings, SimResult,
  PARAMS_SYM, PARAMS_MIS,
  evol, fm, MeV, run_simulation, kadr

# --- SEKCJA 1: STRUKTURY DANYCH (API) ---

"Jednostki"
const fm = 1.0
const MeV = 1 / fm / 197


"""
    HydroParams(C_τπ, C_η, C_λ1)

Niezmienna struktura przechowująca stałe fizyczne modelu hydrodynamiki (BRSSS).
Jej pól nie można zmienić po utworzeniu obiektu.
"""
struct HydroParams
  C_τπ::Float64
  C_η::Float64
  C_λ1::Float64
end

# --- Predefiniowane, niezmienne zestawy parametrów (służą jako "nazwane defaulty") ---
"Parametry z teorii N=4 SYM (holografia/AdS-CFT), pełny model BRSSS."
const PARAMS_SYM = HydroParams(
  (2 - log(2)) / (2 * π),  # C_τπ
  1 / (4 * π),            # C_η
  1 / (2 * π)             # C_λ1
)

"Parametry z pracy PRL (Heller, Spaliński et al.), model zabawkowy MIS."
const PARAMS_MIS = HydroParams(
  (2 - log(2)) / (2 * π),    # C_τπ
  1 / (4 * π),               # C_η
  0.0                        # C_λ1
)


"""
    SimSettings(; n_points=200, τ_start=0.2, ...)

Niezmienna struktura przechowująca ustawienia symulacji.
Posiada konstruktor z argumentami kluczowymi, co pozwala na łatwe
nadpisywanie wartości domyślnych.
"""
struct SimSettings
    ode::Function  
    params::HydroParams
    n_points::Int
    tspan::Tuple{Float64,Float64}
    T_range::Tuple{Float64,Float64}
    A_range::Tuple{Float64,Float64}
end

# Zewnętrzny konstruktor z argumentami kluczowymi. To jest POPRAWNA implementacja
# Twojej intencji "domyślnych ustawień, które można zmieniać".
function SimSettings(;
        ode=ode_brs3!,
        params=PARAMS_SYM,
        n_points=200,
        tspan = (0.2, 10),
        T_range=(300.0 * MeV, 500.0 * MeV),
        A_range=(3.0, 7.0)
    )
  # Wywołuje domyślny, wewnętrzny konstruktor struktury z podanymi wartościami.
  return SimSettings(ode, params, n_points, tspan, T_range, A_range)
end

"""
    SimResult(solutions, settings)

Struktura przechowująca kompletne wyniki symulacji wraz z metadanymi.
"""
struct SimResult
  solutions::Vector{ODESolution}
  settings::SimSettings
end

# --- SEKCJA 2: RDZEŃ SYMULACJI I WIZUALIZACJI ---

# Funkcja wewnętrzna, nieeksportowana
# Definiuje ewolucje w konforemnej teorii BRS3
function ode_brs3!(du, u, p::HydroParams, τ)
  T, A = u
  C_τπ, C_η, C_λ1 = p.C_τπ, p.C_η, p.C_λ1
  du[1] = (T / τ) * (-1 / 3 + A / 18)
  term_T = τ * T * (A + (C_λ1 / (12 * C_η)) * A^2)
  term_A2 = (2 / 9) * C_τπ * A^2
  du[2] = (1 / (C_τπ * τ)) * (8 * C_η - term_T - term_A2)
end

"generuje listę warunków początkowych"
function initial_conditions(settings, seed=5)
    rng = Xoshiro(seed)
    Tmin, Tmax = settings.T_range
    Amin, Amax = settings.A_range
    Ts = rand(rng, Uniform(Tmin, Tmax),settings.n_points)
    As = rand(rng, Uniform(Amin, Amax),settings.n_points)
    return [[Ts[j], As[j]] for j = 1:settings.n_points] 
end


"Ewoluuje dany warunek początkowy"
function evol(u0, p)
  prob = ODEProblem(p.ode, u0, p.tspan, p.params)
  return solve(prob, Tsit5())
end


"""
    run_simulation(params::HydroParams, settings::SimSettings) -> SimResult

Uruchamia pełną symulację i zwraca jeden obiekt z wynikami.
"""
function run_simulation(settings::SimSettings)
  # println("--- Rozpoczynanie Obliczeń Numerycznych...")
  ic = initial_conditions(settings)
  solutions = ODESolution[]
  for u0 in ic
    sol = evol(u0, settings)
    push!(solutions, sol)
  end
  # println("--- Obliczenia Zakończone. ---")
  return SimResult(solutions, settings)
end


"Zwraca (T(t), A(t))"
function TA(solutions, t)
    Ts = [sol(t)[1] for sol in solutions]
    As = [sol(t)[2] for sol in solutions]
    return (Ts, As)
end

function kadr(simres, t)
    (Ts, As) = TA(simres.solutions,t)
    p = plot(title="Kadr w przestrzeni fazowej, t = $t",
             xlabel="Temperatura T [MeV]",
             ylabel="Anizotropia A")
    plot!(p, Ts, As, seriestype=:scatter, label="")
    display(p)
end


end # koniec modułu modHydroSim
