
using Plots
gr()
using LaTeXStrings

# === Przygotowanie danych ===
T_vals = Float64[]
τ_vals = Float64[]
log_ratios = Float64[]  # log(A / A_attr)

for sol in solutions
    for (i, τ) in enumerate(sol.t)
        T_val = sol[1, i]
        A_val = sol[2, i]

        # Pomijamy niefizyczne przypadki (dzielenie przez 0 lub A <= 0)
        if T_val <= 0 || A_val <= 0
            continue
        end

        A_attr = 8*C_η / (τ * T_val) + (16*C_η * C_τπ) / (3 * (τ * T_val)^2)

        # Unika log(0) i dzielenia przez 0
        if A_attr <= 0
            continue
        end

        log_ratio = log(A_val / A_attr)
        push!(T_vals, T_val)
        push!(τ_vals, τ)
        push!(log_ratios, log_ratio)
    end
end

# === Wykres 3D ===
scatter3d(
    T_vals, τ_vals, log_ratios,
    c=τ_vals,
    xlabel="Temperatura T [MeV]",
    ylabel="Czas τ [fm/c]",
    zlabel=L"\ln\left(\frac{A}{A_{\mathrm{attr}}}\right)",
    title="Zbieżność do atraktora hydrodynamicznego (logarytmiczna)",
    markersize=2,
    legend=false
)
