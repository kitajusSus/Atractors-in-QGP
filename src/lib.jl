# Plik: modHydroSim.jl
module modHydroSim

# --- Zależności ---
using DifferentialEquations
using Random

# --- Publiczny Interfejs Modułu ---
export HydroParams, SimSettings, SimResult,
  PARAMS_SYM, PARAMS_MIS,
  evol, fm, MeV, run_simulation

# --- SEKCJA 1: STRUKTURY DANYCH (API) ---

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

"Jednostki"
const fm = 1.0
const MeV = 1 / fm / 197


"""
    SimSettings(; n_points=200, τ_start=0.2, ...)

Niezmienna struktura przechowująca ustawienia symulacji.
Posiada konstruktor z argumentami kluczowymi, co pozwala na łatwe
nadpisywanie wartości domyślnych.
"""
struct SimSettings
  n_points::Int
  τ_start::Float64
  τ_end::Float64
  tspan::Tuple{Float64,Float64}
  T_range::Tuple{Float64,Float64}
  A_range::Tuple{Float64,Float64}
end

# Zewnętrzny konstruktor z argumentami kluczowymi. To jest POPRAWNA implementacja
# Twojej intencji "domyślnych ustawień, które można zmieniać".
function SimSettings(;
  n_points=200,
  τ_start=0.2 * fm,
  τ_end=12 * fm,
  T_range=(300.0 * MeV, 500.0 * MeV),
  A_range=(3.0, 7.0)
)
  # Wywołuje domyślny, wewnętrzny konstruktor struktury z podanymi wartościami.
  return SimSettings(n_points, τ_start, τ_end, (τ_start, τ_end), T_range, A_range)
end

"""
    SimResult(solutions, sol_attractor, params, settings)

Struktura przechowująca kompletne wyniki symulacji wraz z metadanymi.
"""
struct SimResult
  solutions::Vector{ODESolution}
  params::HydroParams
  settings::SimSettings
end

# --- SEKCJA 2: RDZEŃ SYMULACJI I WIZUALIZACJI ---

# Funkcja wewnętrzna, nieeksportowana
function _hydro_evolution!(du, u, p::HydroParams, τ)
  T, A = u
  C_τπ, C_η, C_λ1 = p.C_τπ, p.C_η, p.C_λ1
  du[1] = (T / τ) * (-1 / 3 + A / 18)
  term_T = τ * T * (A + (C_λ1 / (12 * C_η)) * A^2)
  term_A2 = (2 / 9) * C_τπ * A^2
  du[2] = (1 / (C_τπ * τ)) * (8 * C_η - term_T - term_A2)
end

"generuje listę warunków początkowych"
function initial_conditions(settings, seed=5)
  rng = Xoshiro(seed) #seed
  return [
    [rand(rng, settings.T_range[1]:0.1:settings.T_range[2]),
      rand(rng, settings.A_range[1]:0.1:settings.A_range[2])]
    for _ in 1:settings.n_points
  ]
end

"Ewoluuje dany warunek początkowy"
function evol(u0, tspan, params)
  prob = ODEProblem(_hydro_evolution!, u0, tspan, params)
  return solve(prob, Tsit5(), saveat=0.01)
end


"""
    run_simulation(params::HydroParams, settings::SimSettings) -> SimResult

Uruchamia pełną symulację i zwraca jeden obiekt z wynikami.
"""
function run_simulation(params::HydroParams, settings::SimSettings)
  println("--- Rozpoczynanie Obliczeń Numerycznych...")
  solutions = ODESolution[]
  for u0 in initial_conditions
    sol = evol(u0, settings.tspan, params)
    # prob = ODEProblem(_hydro_evolution!, u0, settings.tspan, params)
    # sol = solve(prob, Tsit5(), saveat=0.01)
    push!(solutions, sol)
  end
  println("--- Obliczenia Zakończone. ---")
  return SimResult(solutions, params, settings)
end

end # koniec modułu modHydroSim
