# Krok 1: Importowanie potrzebnych pakietów
using DifferentialEquations
using Plots
using Random
using LaTeXStrings
# Ustawienie backendu dla Plots.jl
gr()

# --- Parametry Fizyczne i Definicja Modelu (bez zmian) ---
C_τπ = 1.0
C_η = 0.75
C_λ1 = 0.0
p = (C_τπ, C_η, C_λ1)

function hydro_evolution!(du, u, p, τ)
    T, A = u
    C_τπ, C_η, C_λ1 = p
    du[1] = (T / τ) * (-1/3 + A / 18)
    term_T = τ * T * (A + (C_λ1 / (12 * C_η)) * A^2)
    term_A2 = (2/9) * C_τπ * A^2
    du[2] = (1 / (C_τπ * τ)) * (8 * C_η - term_T - term_A2)
end

# --- Ustawienia Symulacji (bez zmian) ---
const N_POINTS = 200
τ_start = 0.2
τ_end = 1 #by lepiej zobaczyć reżim wykładniczy, można zmienić na  3.0
tspan = (τ_start, τ_end)
T_min, T_max = 300.0, 500.0
A_min, A_max = 0.0, 4.0
Random.seed!(123)
initial_conditions = [[rand(T_min:0.1:T_max), rand(A_min:0.1:A_max)] for _ in 1:N_POINTS]

# --- Rozwiązywanie Równań ---
println("Rozwiązywanie równań...")
solutions = []
for u0 in initial_conditions
    prob = ODEProblem(hydro_evolution!, u0, tspan, p)
    sol = solve(prob, Tsit5(), saveat=0.01)
    push!(solutions, sol)
end

# --- Wyznaczenie trajektorii atraktora ---
# Startujemy z warunku bliskiego teoretycznemu atraktorowi (regularnemu w τ=0)
# To jest przybliżenie, ale wystarczająco dobre do wizualizacji
u0_attractor = [500.0, 0.1] # T=500 MeV, A=0.1 (blisko regularnego atraktora)
# Dłuższy czas, żeby atraktor się dobrze ustabilizował
prob_attractor = ODEProblem(hydro_evolution!, u0_attractor, (0.01, τ_end), p) 
sol_attractor = solve(prob_attractor, Tsit5(), saveat=0.01)
println("Zakończono rozwiązywanie.")



# === WIZUALIZACJA 1: Poprawiona animacja analityczna ===
T_plot_min = 50
T_plot_max = 550
A_plot_min = -0.5  # Ustawiam limity dla osi Y, by wykres był stabilny
A_plot_max = 2.0

anim = @animate for τ_anim in τ_start:0.005:τ_end
    # Stworzenie pustego wykresu dla każdej klatki
    plot(
        xlabel="Temperatura T [MeV]",
        ylabel=L"ln(\frac{A_{i}}{A_{mis}})", # Poprawiono etykietę dla jasności
        title="Odległość od analitycznego atraktora w τ = $(round(τ_anim, digits=3)) fm/c",
        xlims=(T_plot_min, T_plot_max),
        ylims=(-0.001, :auto), # Ustawione stałe limity
        legend=:bottomleft, # Włączona legenda
        grid=true,
        gridstyle=:dash,
        gridalpha=0.3
    )
    
    # Rysowanie linii referencyjnej (atraktora analitycznego) na y=0
    # To odpowiada sytuacji, gdzie A_numeryczne = A_analityczne, więc log(1) = 0
    hline!([0], 
           color=:magenta, 
           linestyle=:dash, 
           label="Atraktor Analityczny (MIS 2-gi rząd)")
    
    # Pętla po wszystkich rozwiązaniach (kropkach)
    for sol in solutions
        # Używamy interpolacji, aby znaleźć wartość (T,A) w danym czasie
        T_val, A_val = sol(τ_anim)
        
        # Obliczamy wartość atraktora mis dla temperatury DANEJ KROPKI
        # Ta wartość posłuży jako odniesienie A_attr
        A_attr = 8*C_η / (τ_anim * T_val) + (16*C_η * C_τπ) / (3 * (τ_anim * T_val)^2)

        # Rysujemy punkt. Oś Y to logarytm stosunku wartości numerycznej do analitycznej.
        scatter!([T_val], [log(A_val / A_attr)], 
                  marker_z=T_val, 
                  clims=(100, 500), 
                  markerstrokewidth=0, 
                  markersize=4,
                  label="") # Pusty label w pętli
    end
end

gif(anim, "A_T_animacja_log_5fps.gif", fps = 5) # Zmieniono nazwę pliku i fps
gif(anim, "A_T_animacja_log15fps.gif", fps = 15) # Zapisanie animacji do pliku GIF
gif(anim, "A_T_animacja_log_dokładne.gif", fps = 1) # Zapisanie animacji do pliku GIF
println("Animacja została zapisana jako 'A_T_animacja_log.gif'")
# === WIZUALIZACJA 2: Nowa, analityczna animacja log(odległości) ===
println("Generowanie animacji analitycznej...")
anim_log = @animate for τ_anim in τ_start:0.001:τ_end
    # Pobranie wartości na atraktorze w danym czasie
    A_attr = sol_attractor(τ_anim)[2]
    
    plot(
        xlabel="Temperatura T [MeV]",
        ylabel="Odległość od atraktora |A - A_attr|",
        title=" τ = $(round(τ_anim, digits=2)) fm/c",
        legend=false,
        yaxis=:log, # KLUCZOWA ZMIANA: oś Y w skali logarytmicznej!
        #ylims=(1e-4, 10.0) # Zakres dla osi logarytmicznej
    )
    
    # Wykres atraktora jako linii odniesienia (będzie blisko y=0)
    hline!([1e-4], linestyle=:dash, color=:red, label="Atraktor")

    # Rysowanie odległości każdego punktu od atraktora
    for sol in solutions
        T_val, A_val = sol(τ_anim)
        delta_A = abs(A_val - A_attr) + 1e-9 # Dodajemy epsilon, by uniknąć log(0)
        
        # Kolorujemy punkty wg ich początkowej anizotropii, by zobaczyć, czy ma to wpływ
        initial_A = sol.u[1][2]
        scatter!([T_val], [delta_A], marker_z=initial_A, clims=(A_min, A_max), markerstrokewidth=0, markersize=4)
    end
end

gif(anim_log, "hydro_attractor_log_distance.gif", fps = 15)
println("Animacja analityczna została zapisana jako 'hydro_attractor_log_distance.gif'")# Zapisanie animacji do pliku GIF
