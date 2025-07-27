include("lib.jl")
using .modHydroSim
using Plots
gr()


export run_and_plot, plot_trajectories, create_points_animation
"""
    plot_trajectories(wynik_symulacji::SimResult)

Rysuje statyczne trajektorie na podstawie gotowych wyników symulacji.
"""
function plot_trajectories(wynik_symulacji::SimResult)
    println("Generowanie wykresu trajektorii...")
    p = plot(title="Trajektorie w przestrzeni fazowej", xlabel="Temperatura T [MeV]", ylabel="Anizotropia A", legend=false)
    for sol in wynik_symulacji.solutions
        plot!(p, sol, idxs=(1, 2))
    end
    return p
end

"""
    create_points_animation(wynik_symulacji::SimResult; ...)

Tworzy i zapisuje animację GIF pokazującą ewolucję punktów w czasie.
"""
function create_points_animation(
    wynik_symulacji::SimResult;
    output_filename="ewolucja_punktow.gif",
    fps=2,
    frames=200
)
    println("Rozpoczynam generowanie animacji punktów...")

    t_start, t_end = wynik_symulacji.settings.tspan
    tau_values = range(t_start, stop=t_end, length=frames)

    # @animate
    anim = @animate for (idx, tau) in enumerate(tau_values)
        print("\rGenerowanie klatki $(idx)/$(frames)...") 

        # Stwórz bazowy wykres dla pojedynczej klatki
        p = plot(
            title="Ewolucja punktów (τ = $(round(tau, digits=3)))",
            xlabel="Temperatura T [MeV]",
            ylabel="Anizotropia A",
            # Stałe granice osi są kluczowe dla dobrej animacji!
            xlims=(0.0, wynik_symulacji.settings.T_range[2]),
            ylims=(-1.0, wynik_symulacji.settings.A_range[2] + 1)
        )

        # Dla każdego rozwiązania,
        for (i, sol)in enumerate(wynik_symulacji.solutions)
            aktualny_stan = sol(tau)
            scatter!(p, [aktualny_stan[1]], [aktualny_stan[2]], label="", seriescolor=i, markersize=3)
        end
    end

    # Zapisz animację do pliku GIF
    gif(anim, output_filename, fps=fps)
    println("\nAnimacja została zapisana jako: $(output_filename)")
end


"""
    run_and_plot(settings::SimSettings = SimSettings())

Główna funkcja orkiestrująca: uruchamia symulację i rysuje statyczny wykres.
"""
function run_and_plot(settings::SimSettings = SimSettings())
    println("--- Rozpoczynam symulację... ---")
    wynik = run_simulation(settings)
    println("--- Symulacja zakończona. ---")
    
    wykres = plot_trajectories(wynik)
    display(wykres)
    println("\nWykres statyczny gotowy.")
    
    return (wynik, wykres)
end

println("Plik `test_1.jl` wczytany. Dostępne funkcje: `run_and_plot`, `plot_trajectories`, `create_points_animation`.")
