
using DifferentialEquations
using Plots

# --- Parametry fizyczne ---
const Cη = 1 / (4π)             # eta/s
const Cτ = (2 - log(2)) / (2π)  # tau_pi * T

# --- Równanie różniczkowe dla A(w) ---
function attractor_eq!(dA, A, p, w)
    numerator = (3/2) * (8 * Cη / w - A[1]) - (Cτ / (3w)) * A[1]^2
    denominator = Cτ * (1 + A[1]/12)
    dA[1] = numerator / denominator
end

# --- Zakres zmiennej w (tau * T) ---
wspan = (0.2, 5)

# --- Warunki początkowe (kilka trajektorii) ---
A0_list = [-0.5, 0.0, 1.0, 3.0, 4.0, -50, 10]
sols = []

for A0 in A0_list
    prob = ODEProblem(attractor_eq!, [A0], wspan)
    sol = solve(prob, Tsit5(), abstol=1e-8, reltol=1e-8)
    push!(sols, sol)
end

# --- Tworzenie wykresu ---
plt = plot(size=(800,600))
colors = [:blue, :blue, :blue, :blue, :blue, :blue, :orange]

for (i, sol) in enumerate(sols)
    plot!(plt, sol.t, [u[1] for u in sol.u], label="A₀ = $(A0_list[i])", color=colors[i], lw=2)
end

xlabel!(plt, "w = τ·T")
ylabel!(plt, "A(w)")
title!(plt, "Ewolucja A(w) dla różnych warunków początkowych")

# --- Zapis do pliku PNG ---
savefig(plt, "attractor_plot1.png")
println("Wykres zapisany jako attractor_plot.png")

