
using DifferentialEquations
using Plots
using LaTeXStrings
using Plots
plotlyjs()


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
wspan = (1, 50)

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
colors = [:blue, :green, :blue, :blue, :blue, :red, :blue]

for (i, sol) in enumerate(sols)
    plot!(plt, sol.t, [u[1] for u in sol.u], label="A₀ = $(A0_list[i])", color=colors[i], lw=2)
end

a0_list = [-50 + 0.5*i for i in 1:500]

sols2 = []
for a0 in a0_list
    prob2 = ODEProblem(attractor_eq!,[a0],wspan)
    sol2 = solve(prob2,Tsit5(), abstol=1e-8, reltol = 1e-8)
    push!(sols2,sol2)
end


kolory = [:violet]
# inne warianty
for (i,sol2) in enumerate(sols2)
    plot!(plt, sol2.t, [u[1] for u in sol2.u], label ="a_0 = $(a0_list[i])", color =:violet, lw=2)
    
end

xlabel!(plt, "w = τ·T")
ylabel!(plt, "A(w)")
title!(plt, "Ewolucja A(w) dla różnych warunków początkowych")
display(plt)
# --- Zapis do pliku PNG ---
# savefig(plt, "attractor_plot1.png")
# println("Wykres zapisany jako attractor_plot.png")

