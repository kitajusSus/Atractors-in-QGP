# --- Załadowanie potrzebnych pakietów ---
using DifferentialEquations
using Plots
using LaTeXStrings

# Ustawienie backendu dla interaktywnych wykresów 3D
plotlyjs() 

# --- Parametry fizyczne ---
const Cη = 1 / (4π)             # eta/s
const Cτ = (2 - log(2)) / (2π)  # tau_pi * T

# --- Równanie różniczkowe dla [A(w), T(w)] ---
# u[1] = A(w), u[2] = T(w)
# p to parametry (nieużywane tutaj), w to zmienna niezależna
function attractor_system!(du, u, p, w)
    A = u[1]
    T = u[2]
    
    # 1. Równanie na A'(w)
    # Dodajemy zabezpieczenie, aby uniknąć dzielenia przez zero, gdy A jest bliskie -12
    denominator_A = Cτ * (1 + A/12)
    if abs(denominator_A) < 1e-9
        denominator_A = 1e-9 * sign(denominator_A)
    end
    numerator_A = (3/2) * (8 * Cη / w - A) - (Cτ / (3w)) * A^2
    du[1] = numerator_A / denominator_A

    # 2. Równanie na T'(w)
    # Podobnie zabezpieczamy się przed A ≈ -12
    denominator_T = w * (A + 12)
    if abs(denominator_T) < 1e-9
         denominator_T = 1e-9 * sign(denominator_T)
    end
    du[2] = T * (A - 6) / denominator_T
end

# --- Ustawienia symulacji ---
w₀ = 0.05
w_final = 3.0
wspan = (w₀, w_final)

# Warunki początkowe A₀ w w₀
initial_A = -2.0:1.0:8.0  # Zestaw różnych A₀

# Warunek początkowy T₀ w w₀ (wspólny dla wszystkich trajektorii)
T₀ = 1.0 # Możesz go zmieniać, np. na 0.5, 2.0 etc.

# --- Rozwiązywanie układu równań ---
println("Rozwiązywanie układu równań dla $(length(initial_A)) warunków początkowych...")
sols = []
for A₀ in initial_A
    u₀ = [A₀, T₀] # Wektor stanu początkowego
    prob = ODEProblem(attractor_system!, u₀, wspan)
    sol = solve(prob, Tsit5(), abstol=1e-8, reltol=1e-8, saveat=0.01)
    push!(sols, sol)
end
println("Rozwiązywanie zakończone.")


# --- Ekstrakcja danych i obliczanie τ ---
# Przygotowujemy dane do wykresu 3D
w_vals = sols[1].t
trajectories_3d = []

for sol in sols
    # sol.u to wektor wektorów [[A1, T1], [A2, T2], ...]
    A_vals = [u[1] for u in sol.u]
    T_vals = [u[2] for u in sol.u]
    
    # Obliczamy τ(w) = w / T(w)
    τ_vals = sol.t ./ T_vals
    
    # Przechowujemy trójki (w, τ, A)
    push!(trajectories_3d, (sol.t, τ_vals, A_vals))
end


# --- Tworzenie statycznego wykresu 2D (A vs w) ---
# Używamy backendu GR do szybkiego zapisu
plotlyjs()
plt_2d = plot!(
    size=(800, 600),
    title="Ewolucja A(w) w kierunku atraktora (2D)",
    xlabel=L"w = \tau T",
    ylabel=L"A(w)",
    legend=:topright,
    ylims=(-5, 10),
    xlims=(w₀, w_final)
)

# Rysowanie wszystkich trajektorii A(w)
for sol in sols
    plot!(plt_2d, sol, vars=(0,1), label="", color=:blue, alpha=0.6)
end

# Znajdujemy i wyróżniamy atraktor
# Atraktorem jest rozwiązanie, które startuje z odpowiedniej wartości A₀
# Możemy go przybliżyć jako rozwiązanie, które jest najniżej dla małych w
# W tym przypadku, to często ostatnie rozwiązanie na naszej liście startującej od ujemnych A
# Lub możemy go zdefiniować przez specjalny warunek początkowy, jeśli go znamy.
# Tutaj, dla prostoty, wyróżnimy jedną z trajektorii jako "reprezentanta" atraktora.
idx_attractor = 1 # Pierwsze rozwiązanie (startujące z A₀ = -2.0)
plot!(plt_2d, sols[idx_attractor], vars=(0,1), color=:red, lw=3, label="Trajektoria atraktora")

#display(plt_2d)
#savefig(plt_2d, "attractor_2d_plot.png")
println("Zapisano wykres 2D do pliku attractor_2d_plot.png")


# --- Tworzenie interaktywnego wykresu 3D (A vs w vs τ) ---
# Wracamy do backendu plotlyjs() dla interaktywności
plotlyjs()
plt_3d = plot!(
    size=(1000, 800),
    title="Ewolucja w przestrzeni fazowej (w, τ, A)",
    xlabel=L"w = \tau T",
    ylabel=L"\tau \text{ (czas własny)}",
    zlabel=L"A(w) \text{ (anizotropia)}",
    camera=(30, 30), # Ustawienie początkowe kamery (kąt, elewacja)
    legend=false
)

# Rysowanie trajektorii w 3D
for traj in trajectories_3d
    # traj to krotka (w_vals, τ_vals, A_vals)
    plot!(plt_3d, traj[1], traj[2], traj[3], color=:blue, alpha=0.7)
end

# Wyróżnienie trajektorii atraktora w 3D
attractor_traj_3d = trajectories_3d[idx_attractor]
plot!(plt_3d, attractor_traj_3d[1], attractor_traj_3d[2], attractor_traj_3d[3], color=:red, lw=3)

display(plt_3d)
println("Wygenerowano interaktywny wykres 3D. Możesz go obracać i przybliżać.")
