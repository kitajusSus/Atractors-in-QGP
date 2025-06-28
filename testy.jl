# --- Załadowanie potrzebnych pakietów ---
using DifferentialEquations
using Plots
using LaTeXStrings

# --- Ustawienie backendu dla Plots.jl ---
# gr() jest znacznie szybszy do generowania animacji i zapisu plików.
gr() 

# --- Parametry fizyczne ---
const Cη = 1 / (4π)             # eta/s
const Cτ = (2 - log(2)) / (2π)  # tau_pi * T

# --- Równanie różniczkowe dla A(w) ---
# Definiuje dynamikę systemu. dA to wektor pochodnych, A to wektor stanu.
function attractor_eq!(dA, A, p, w)
    # Równanie ma osobliwość, gdy mianownik jest bliski zera.
    # A[1] ≈ -12. Solver powinien sobie z tym poradzić.
    numerator = (3/2) * (8 * Cη / w - A[1]) - (Cτ / (3w)) * A[1]^2
    denominator = Cτ * (1 + A[1]/12)
    dA[1] = numerator / denominator
end

# --- Zakres zmiennej w (tau * T) ---
wspan = (0.1, 2)

# --- Warunki początkowe ---
# Wybrane trajektorie do wyróżnienia
A0_highlight = [-50.0, -0.5, 0.0, 1.0, 4.0, 10.0]
# Trajektorie tła, pokazujące ogólny przepływ
A0_background = range(-100, 20, length=100)
# Połączenie wszystkich warunków początkowych
all_A0 = vcat(A0_highlight, A0_background)

# --- Rozwiązywanie równań różniczkowych ---
println("Rozwiązywanie równań dla $(length(all_A0)) warunków początkowych...")
sols = []
for A0 in all_A0
    prob = ODEProblem(attractor_eq!, [A0], wspan)
    sol = solve(prob, Tsit5(), abstol=1e-8, reltol=1e-8, saveat=0.05)
    push!(sols, sol)
end
println("Rozwiązywanie zakończone.")

# --- Tworzenie statycznego wykresu (dla wglądu i zapisu) ---
println("Tworzenie statycznego wykresu...")
# Ustawienie motywu i rozmiaru
plt_static = plot(
    size=(800, 600),
    theme=:solarized,
    title="Ewolucja A(w) w kierunku atraktora",
    xlabel=L"w = \tau \cdot T",
    ylabel=L"A(w)",
    legend=:topright,
    ylims=(-15, 15) # Stałe granice dla spójności
)

# Kolory dla wyróżnionych trajektorii
highlight_colors = [:red, :blue, :green, :orange, :purple, :cyan]
num_highlight = length(A0_highlight)

# Rysowanie wyróżnionych trajektorii
for i in 1:num_highlight
    plot!(plt_static, sols[i],
        label=L"A_0 = %$(A0_highlight[i])",
        color=highlight_colors[i],
        lw=2.5 # Grubsza linia
    )
end

# Rysowanie trajektorii tła
for i in (num_highlight + 1):length(sols)
    plot!(plt_static, sols[i],
        label="", # Brak etykiety w legendzie
        color=:black,
        alpha=0.3, # Przezroczystość
        lw=1 # Cieńsza linia
    )
end

# Wyświetlenie wykresu statycznego
display(plt_static)

# Zapis do pliku PNG
savefig(plt_static, "attractor_plot_static.png")
println("Statyczny wykres zapisany jako attractor_plot_static.png")


# --- Tworzenie animacji ---
println("Tworzenie animacji (może to chwilę potrwać)...")

# `anim` to obiekt, do którego będziemy dodawać klatki
anim = @animate for w_current in range(wspan[1] + 0.01, wspan[2], length=200)
    
    # Tworzymy klatkę - czysty wykres z tytułem pokazującym postęp
    p_anim = plot(
        size=(800, 600),
        theme=:solarized,
        title=L"Ewolucja do $w = %$(round(w_current, digits=2))$",
        xlabel=L"w = \tau \cdot T",
        ylabel=L"A(w)",
        legend=:topright,
        #ylims=(-20, 15) # Kluczowe dla stabilnej animacji!
    )

    # Rysowanie wyróżnionych trajektorii do punktu w_current
    for i in 1:num_highlight
        # `tspan` ogranicza rysowanie rozwiązania do danego przedziału
        plot!(p_anim, sols[i],
            tspan=(wspan[1], w_current),
            label=L"A_0 = %$(A0_highlight[i])",
            color=highlight_colors[i],
            lw=2.5
        )
    end
    
    # Rysowanie trajektorii tła do punktu w_current
    for i in (num_highlight + 1):length(sols)
        plot!(p_anim, sols[i],
            tspan=(wspan[1], w_current),
            label="",
            color=:black,
            alpha=0.3,
            lw=1
        )
    end
end

# Zapis animacji do pliku GIF
gif(anim, "attractor_animation.gif", fps=30)
println("Animacja zapisana jako attractor_animation.gif")
