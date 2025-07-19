# --- Załadowanie potrzebnych pakietów ---
using DifferentialEquations
using Plots
using LaTeXStrings
using Colors, ColorSchemes # Dodajemy pakiety do obsługi kolorów

# Ustawienie backendu dla interaktywnych wykresów 3D
plotlyjs()

# --- Parametry fizyczne ---
const Cη = 1 / (4π)
const Cτ = (2 - log(2)) / (2π)
const ħc = 197.3 # MeV * fm

# --- Równanie różniczkowe dla [A(w), T(w)] ---
function attractor_system!(du, u, p, w)
    A, T = u

    denominator_A = Cτ * (1 + A / 12)
    denominator_A = abs(denominator_A) < 1e-9 ? 1e-9 * sign(denominator_A) : denominator_A

    numerator_A = (3/2) * (8 * Cη / w - A) - (Cτ / (3w)) * A^2
    du[1] = numerator_A / denominator_A

    denominator_T = w * (A + 12)
    denominator_T = abs(denominator_T) < 1e-9 ? 1e-9 * sign(denominator_T) : denominator_T
    du[2] = T * (A - 6) / denominator_T
end

# --- Ustawienia symulacji ---
w₀_start = 0.05
w_final = 3.0
wspan = (w₀_start, w_final)

initial_A = -11.5:1.0:12.0
T₀ = 1.0

# --- Rozwiązywanie układu równań ---
println("Rozwiązywanie układu równań dla $(length(initial_A)) warunków początkowych...")
sols = []
for A₀ in initial_A
    # Przeliczamy T₀ na odpowiednią wartość w₀
    # Zauważ, że T₀ jest w jednostkach 1/fm, a nie MeV tutaj
    # Jeśli chcesz zadać T₀ w MeV, musisz to uwzględnić. Załóżmy, że T₀ jest już w 1/fm.
    w₀ = w₀_start # Startujemy wszystkie trajektorie z tego samego w
    u₀ = [A₀, T₀ * (w₀/w₀_start)] # Skalujemy T₀, by odpowiadało w₀
    
    prob = ODEProblem(attractor_system!, u₀, wspan)
    # Rodas5 jest dobrym wyborem dla problemów sztywnych (stiff)
    sol = solve(prob, Rodas5(), abstol=1e-8, reltol=1e-8, saveat=0.01)
    push!(sols, sol)
end
println("Rozwiązywanie zakończone.")

# --- Ekstrakcja danych i obliczanie τ ---
trajectories_3d = []
for sol in sols
    A_vals = [u[1] for u in sol.u]
    T_vals = [u[2] for u in sol.u]
    # Teraz T jest w jednostkach 1/fm, więc nie trzeba mnożyć przez ħc
    τ_vals = sol.t ./ T_vals 
    push!(trajectories_3d, (sol.t, τ_vals, A_vals))
end

# --- Tworzenie wykresu 2D (A vs w) ---
println("Tworzenie wykresu 2D...")
gr() # Przełączamy na GR dla szybkiego zapisu statycznego
plt_2d = plot(
    size=(800, 600),
    title="Ewolucja A(w) w kierunku atraktora (2D)",
    xlabel=L"w = \tau T",
    ylabel=L"A(w)",
    legend=:topright,
    ylims=(-5, 10),
    xlims=wspan,
    grid=true
)

# Generowanie gradientu kolorów
gradient_kolorow = cgrad(:viridis, length(sols), categorical=true)

# Rysowanie trajektorii
for (i, sol) in enumerate(sols)
    plot!(plt_2d, sol.t, [u[1] for u in sol.u], label="", color=gradient_kolorow[i], alpha=0.7, lw=2)
end

# Wyróżnienie trajektorii atraktora
# Wybór środkowej trajektorii jest rozsądnym przybliżeniem
id_attractor = findfirst(A -> A ≈ 0.0, initial_A)
if isnothing(id_attractor)
    id_attractor = Int(cld(length(sols), 2))
end
plot!(plt_2d, sols[id_attractor].t, [u[1] for u in sols[id_attractor].u], color=:red, lw=3, label="Atraktor")

display(plt_2d)
savefig(plt_2d, "attractor_2d_plot.png")
println("Zapisano wykres 2D do pliku attractor_2d_plot.png")

# --- Tworzenie wykresu 3D (A vs w vs τ) ---
println("Tworzenie wykresu 3D...")
plotlyjs() # Wracamy do interaktywnego backendu
plt_3d = plot(
    size=(1000, 800),
    title="Ewolucja w przestrzeni fazowej (w, τ, A)",
    xlabel="w = τT",
    ylabel="τ (czas własny)",
    zlabel="A(w) (anizotropia)",
    camera=(30, 30),
    legend=false,
    grid=true
)

# Rysowanie wszystkich trajektorii w 3D
for (i, traj) in enumerate(trajectories_3d)
    plot!(plt_3d, traj[1], traj[2], traj[3], color=gradient_kolorow[i], alpha=0.8, lw=4)
end

# Wyróżnienie atraktora
attractor_traj_3d = trajectories_3d[id_attractor]
plot!(plt_3d, attractor_traj_3d[1], attractor_traj_3d[2], attractor_traj_3d[3], color=:red, lw=5)

display(plt_3d)
println("Wygenerowano interaktywny wykres 3D.")
