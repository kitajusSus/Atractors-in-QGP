include("lib.jl")
using .modHydroSim
using Plots
using MultivariateStats

"""
    prepare_hydro_data(sim_settings::SimSettings, n_time_samples::Int) -> Matrix, Vector
przekształca trajektorie w dane wysokowymiarowe gotowe do użycia przez PCA.
"""
function prepare_hydro_data(sim_settings::SimSettings, n_time_samples::Int)
    println("--- Krok 1: Generowanie danych z symulacji hydrodynamiki... ---")
    wynik_symulacji = run_simulation(sim_settings)
    n_trajectories = length(wynik_symulacji.solutions)
    
    t_start, t_end = wynik_symulacji.settings.tspan
    sample_times = range(t_start, stop=t_end, length=n_time_samples)
    D = n_time_samples * 2
    X = zeros(Float64, n_trajectories, D)
    initial_temperatures = zeros(Float64, n_trajectories)

    println("--- Krok 2: Przekształcanie trajektorii w wektory wysokowymiarowe... ---")
    for i in 1:n_trajectories
        sol = wynik_symulacji.solutions[i]
        sampled_points = sol(sample_times)
        flat_vector = vcat([[p[1], p[2]] for p in sampled_points]...)
        X[i, :] = flat_vector
        initial_temperatures[i] = sol.u[1][1]
    end
    
    println("Przygotowanie danych zakończone. Wymiar danych wejściowych: $n_trajectories x $D")
    return X, initial_temperatures
end

function main()
    hydro_settings = SimSettings(
        n_points=250, 
        tspan=(0.2, 0.8)
    )
    n_time_samples = 5 # Liczba próbek z każdej trajektorii

    X, initial_temps = prepare_hydro_data(hydro_settings, n_time_samples)

    println("\n--- Krok 3: Redukcja wymiarowości PCA ---")
    pca_model = fit(PCA, X; maxoutdim=2)    # 2-wymiarowa projekcja
    Y = transform(pca_model, X)

    println("\n--- Krok 4: Generowanie wykresu wynikowego... ---")
    p = scatter(
        Y[:, 1], Y[:, 2],
        marker_z=initial_temps,
        title="Wizualizacja PCA trajektorii hydrodynamiki",
        xlabel="Komponent PCA 1",
        ylabel="Komponent PCA 2",
        label="",
        markersize=5,
        alpha=0.8,
        color=:plasma,
        colorbar_title="  Temp. początkowa [MeV]"
    )
    display(p)
    println("Gotowe. Wyświetlono wykres.")
end

main()
