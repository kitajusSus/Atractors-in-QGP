using AttractorsQGP
using GLMakie
using Statistics
using LinearAlgebra

function run_hjsw_lpca_example()
    println("=== Running HJSW Model LPCA and Attractor Analysis ===")

    # 1. Initialize HJSW Model (N=4 SYM parameters)
    model = HJSWModel(
        eta_over_s = 1 / (4 * π),
        omega_R = 9.8,
        omega_I = 8.629
    )

    # 2. Generate initial conditions in 3D: [T0, A0, B0]
    n_trajectories = 50
    ics = generate_initial_conditions(
        model,
        n_trajectories;
        T_range = (400.0, 1200.0),
        A_range = (-2.0, 8.0),
        B_range = (-4.0, 4.0),
        seed = 42
    )

    println("Generated $(length(ics)) 3D initial conditions [T, A, B].")

    # 3. Solve Bjorken hydrodynamics trajectories
    tspan = (0.2, 1.5)
    sols = generate_trajectories(model, ics, tspan; saveat = 0.01, parallel = :serial)
    dataset = build_dataset(sols; temperature_unit = :fm)
    println("Dataset built with shape: $(size(dataset)) (columns: [tau, T, A, B])")

    # 4. Perform LPCA analysis on feature columns [T, A, B]
    n_slices = 15
    tablica_k = [5, 10, 15, 20]

    println("Running Local PCA analysis across time slices...")
    fig_lpca = plot_local_pca(
        dataset;
        n_slices = n_slices,
        tablica_k = tablica_k,
        feature_cols = [2, 3, 4]
    )

    # 5. Plot 3D Phase Space Evolution (T, A, B)
    println("Plotting 3D Phase Space Evolution...")
    fig_phase_space_3d = plot_phase_space_evolution_3d(dataset)

    # 6. Plot Attractor (w = tau*T vs A)
    println("Plotting 2D Attractor...")
    fig_attractor = plot_attractor(dataset; n_trajectories = 15)

    # 7. Save generated plots
    output_dir = joinpath(@__DIR__, "..", "..", "plots")
    mkpath(output_dir)

    save(joinpath(output_dir, "hjsw_lpca_dimension.png"), fig_lpca)
    save(joinpath(output_dir, "hjsw_phase_space_3d.png"), fig_phase_space_3d)
    save(joinpath(output_dir, "hjsw_attractor.png"), fig_attractor)

    println("All plots successfully saved to: ", output_dir)
    return (fig_lpca = fig_lpca, fig_phase_space_3d = fig_phase_space_3d, fig_attractor = fig_attractor)
end

if abspath(PROGRAM_FILE) == @__FILE__
    run_hjsw_lpca_example()
end
