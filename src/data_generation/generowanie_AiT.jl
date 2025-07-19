# ===================================================================
#     GENERATOR DANYCH QGP - WERSJA "SLOW-MOTION START"
#
#   Ten skrypt generuje dane z bardzo wysoką rozdzielczością
#   na początku, aby płynnie zobrazować gwałtowny spadek anizotropii.
# ===================================================================

using DifferentialEquations, CSV, DataFrames, Printf

println("Rozpoczynam generowanie danych z ultra-wysoką rozdzielczością na starcie.")

# Definicje fizyczne (bez zmian)
const C_η = 1 / (4π); const C_τπ = (2 - log(2)) / (2π); const C_λ1 = 1 / (2π)
const PARAMS = (C_η, C_τπ, C_λ1)

function hydro_evolution!(du, u, p, τ)
    T, A = u
    C_η, C_τπ, C_λ1 = p
    du[1] = (T / τ) * (-1/3 + A / 18)
    term_T_dep = τ * T * (A + (C_λ1 / (12 * C_η)) * A^2)
    term_A2_dep = (2/9) * C_τπ * A^2
    du[2] = (1 / (C_τπ * τ)) * (8*C_η - term_T_dep - term_A2_dep)
end

# Parametry symulacji
const N_POINTS = 500
const T_RANGE = (200, 400) # Wyższe T dla dramatyzmu
const A_RANGE = (2, 5)   # Wysokie A, żeby dobrze widzieć spadek

const τ_start = 0.2
const τ_end = 1.2
const τ_span = (τ_start, τ_end)

# --- KLUCZOWA ZMIANA: ADAPTACYJNE PRÓBKOWANIE CZASU ---
# Chcemy zobaczyć spadek A w "zwolnionym tempie".
# Krok 1: Bardzo gęsto na samym początku (pierwsze 0.05 fm/c)
save_points_early = range(τ_start, 0.25, length=50)
# Krok 2: Mniej gęsto w dalszej części
save_points_late = range(0.25 + 0.01, τ_end, length=51)
# Łączymy je w jedną listę punktów do zapisu
const τ_save_points = unique(sort([save_points_early..., save_points_late...]))
# --------------------------------------------------------

println("Zastosowano adaptacyjne próbkowanie. Pierwsze 5 kroków czasowych:")
for i in 1:5
    @printf "τ_%d = %.5f\n" i τ_save_points[i]
end
println("...")

# Inicjalizacja DataFrame (bez zmian)
df = DataFrame(
    Run_ID = Int[], T_0 = Float64[], A_0 = Float64[],
    tau = Float64[], T_at_tau = Float64[], A_at_tau = Float64[]
)

# Główna pętla symulacji
println("Przeprowadzanie symulacji...")
initial_conditions_T = T_RANGE[1] .+ (T_RANGE[2] - T_RANGE[1]) .* rand(N_POINTS)
initial_conditions_A = A_RANGE[1] .+ (A_RANGE[2] - A_RANGE[1]) .* rand(N_POINTS)

for i in 1:N_POINTS
    if i % 25 == 0; @printf "Symulacja %d / %d\n" i N_POINTS; end
    u₀ = [initial_conditions_T[i], initial_conditions_A[i]]
    prob = ODEProblem(hydro_evolution!, u₀, τ_span, PARAMS)
    sol = solve(prob, Tsit5(), saveat=τ_save_points)
    for (j, τ) in enumerate(sol.t)
        push!(df, (i, u₀[1], u₀[2], τ, sol.u[j][1], sol.u[j][2]))
    end
end

# Zapisanie do pliku CSV
output_filename = "wyniki_symulacji_QGP2.csv"
println("\nZapisywanie wyników do pliku: $output_filename")
CSV.write(output_filename, df)
println("\nGotowe!")
