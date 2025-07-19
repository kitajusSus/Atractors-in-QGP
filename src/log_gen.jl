# Logarytmiczne próbkowanie czasowe w symulacji hydrodyamicznej
# 1. Importowanie potrzebnych pakietów
using DifferentialEquations
using CSV
using DataFrames
using Printf

println("Rozpoczynam generowanie danych z logarytmicznym próbkowaniem w czasie.")

# 2. Definicja stałych fizycznych i modelu (bez zmian)
const C_η = 1 / (4π)
const C_τπ = (2 - log(2)) / (2π)
const C_λ1 = 1 / (2π)
const PARAMS = (C_η, C_τπ, C_λ1)

function hydro_evolution!(du, u, p, τ)
    T, A = u
    C_η, C_τπ, C_λ1 = p
    du[1] = (T / τ) * (-1/3 + A / 18)
    term_T_dep = τ * T * (A + (C_λ1 / (12 * C_η)) * A^2)
    term_A2_dep = (2/9) * C_τπ * A^2
    du[2] = (1 / (C_τπ * τ)) * (8*C_η - term_T_dep - term_A2_dep)
end

# 3. Parametry symulacji
const N_POINTS = 500
const T_RANGE = (200.0, 500.0)
const A_RANGE = (-1.0, 5.0)

const τ_start = 0.2
const τ_end = 1.2
const τ_span = (τ_start, τ_end)

# --- KLUCZOWA ZMIANA: LOGARYTMICZNE PRÓBKOWANIE CZASU ---
const N_SAVE_POINTS = 101 # Wciąż chcemy 101 punktów na trajektorię
# Tworzymy punkty na skali logarytmicznej, a potem wracamy do skali liniowej
log_tau_start = log(τ_start)
log_tau_end = log(τ_end)
log_tau_points = range(log_tau_start, log_tau_end, length=N_SAVE_POINTS)
# Konwertujemy z powrotem na punkty w czasie τ. Będą one gęstsze na początku.
const τ_save_points = exp.(log_tau_points)
# --------------------------------------------------------

# Wyświetlmy pierwsze kilka punktów, żeby zobaczyć różnicę:
println("Przykładowe punkty czasowe do zapisu:")
for i in 1:5
    @printf "τ_%d = %.4f\n" i τ_save_points[i]
end
println("...")

# 4. Inicjalizacja DataFrame (bez zmian)
df = DataFrame(
    Run_ID = Int[],
    T_0 = Float64[],
    A_0 = Float64[],
    tau = Float64[],
    T_at_tau = Float64[],
    A_at_tau = Float64[]
)

# 5. Główna pętla symulacji (bez zmian, ale teraz używa τ_save_points)
println("Przeprowadzanie symulacji...")

initial_conditions_T = T_RANGE[1] .+ (T_RANGE[2] - T_RANGE[1]) .* rand(N_POINTS)
initial_conditions_A = A_RANGE[1] .+ (A_RANGE[2] - A_RANGE[1]) .* rand(N_POINTS)

for i in 1:N_POINTS
    if i % 25 == 0; @printf "Symulacja %d / %d\n" i N_POINTS; end

    u₀ = [initial_conditions_T[i], initial_conditions_A[i]]
    prob = ODEProblem(hydro_evolution!, u₀, τ_span, PARAMS)
    
    # Solver zapisze wyniki DOKŁADNIE w naszych logarytmicznie rozmieszczonych punktach
    sol = solve(prob, Tsit5(), saveat=τ_save_points)

    for (j, τ) in enumerate(sol.t)
        T_val, A_val = sol.u[j]
        push!(df, (i, u₀[1], u₀[2], τ, T_val, A_val))
    end
end

# 6. Zapisanie do pliku CSV
output_file = "dane_log_tau.csv"
println("\nZapisywanie wyników do pliku: $output_file")
CSV.write(output_file, df)

println("\nGotowe! Plik '$output_file' został pomyślnie utworzony.")
