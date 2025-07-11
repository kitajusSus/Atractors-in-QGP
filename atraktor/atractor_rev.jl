
using DifferentialEquations
using PlotlyJS


# --- Parametry fizyczne ---
const Cη = 1 / (4π)
const Cτ = (2 - log(2)) / (2π)

# --- Równanie różniczkowe dla A(w) ---
function attractor_eq!(dA, A, p, w)
    numerator = (3/2) * (8 * Cη / w - A[1]) - (Cτ / (3w)) * A[1]^2
    denominator = Cτ * (1 + A[1]/12)
    dA[1] = numerator / denominator
end

# --- Zakres zmiennej w (tau * T), od nieskończoności do zera ---
wspan = (1e4, 0.5)  # od "nieskończoności" (duże w) do małego w

# --- Warunki początkowe (kilka trajektorii) ---
A0_list = [-0.5, 0.0, 1.0, 3.0, 4.0, -50.0, 10.0]
sols = []

for A0 in A0_list
    prob = ODEProblem(attractor_eq!, [A0], wspan)
    sol = solve(prob, Tsit5(), abstol=1e-8, reltol=1e-8)
    push!(sols, sol)
end

# --- Tworzenie śladów dla pierwszej grupy trajektorii ---
traces = GenericTrace[]
colors_A0 = Dict(-0.5 => "blue", 0.0 => "green", 1.0 => "blue", 3.0 => "blue", 4.0 => "blue", -50.0 => "red", 10.0 => "blue")

for (i, sol) in enumerate(sols)
    A0_val = A0_list[i]
    push!(traces, scatter(x=sol.t, y=[u[1] for u in sol.u],
                          mode="lines", name="A₀ = $(A0_val)",
                          line=attr(color=get(colors_A0, A0_val, "blue"), width=2)))
end

# --- Dodatkowe warunki początkowe dla drugiej grupy trajektorii ---
a0_list = [-25 + 2*i for i in 1:25]
sols2 = []

for a0 in a0_list
    prob2 = ODEProblem(attractor_eq!, [a0], wspan)
    sol2 = solve(prob2, Tsit5(), abstol=1e-8, reltol=1e-8)
    push!(sols2, sol2)
end

# --- Tworzenie śladów dla drugiej grupy trajektorii ---
for (i, sol2) in enumerate(sols2)
    push!(traces, scatter(x=sol2.t, y=[u[1] for u in sol2.u],
                          mode="lines", name="a₀ = $(a0_list[i])",
                          line=attr(color="violet", width=2), showlegend=false))
end

# --- Tworzenie układu wykresu ---
my_layout = Layout(
    title="Ewolucja A(w) od nieskończoności do 0",
    xaxis=attr(title="w = τ·T (maleje)", autorange="reversed"),
    yaxis=attr(title="A(w)"),
    height=600,
    width=800
)

# --- Tworzenie i wyświetlanie wykresu ---
plt = plot(traces, my_layout)
display(plt)

# --- Zapis do pliku PNG (opcjonalnie) ---
# savefig(plt, "attractor_plot_backward.png")
# println("Wykres zapisany jako attractor_plot_backward.png")
