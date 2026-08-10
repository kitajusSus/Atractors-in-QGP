using AttractorsQGP
using GLMakie
using Statistics
using LinearAlgebra

"""
    run_hjsw_lpca_regularization_example()

Demonstrates Local PCA (LPCA) analysis on HJSW model state space data (tau, T, A, B)
comparing four different data regularization/normalization methods:
1. Abs-Max scaling (`:max`) -> scales columns to [-1, 1]
2. Min-Max scaling (`:minmax`) -> scales columns to [0, 1]
3. Z-score Standardization (`:zscore`) -> zero mean, unit variance
4. Raw unscaled data (`:none`) -> shows effect of raw physical units
"""
function run_hjsw_lpca_regularization_example()
    println("=== Running HJSW Model LPCA with Data Regularization Comparisons ===")

    # 1. Initialize HJSW Model
    model = HJSWModel(
        eta_over_s = 1 / (4 * π),
        omega_R = 9.80005,
        omega_I = 2.87631
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

    # 4. Generate multi-panel comparison plot across 4 regularization methods
    n_slices = 15
    tablica_k = [5, 10, 15, 20]
    methods = [:max, :minmax, :zscore, :none]

    println("Generating multi-panel LPCA comparison across regularization methods...")
    fig_comparison = plot_local_pca_regularizations(
        dataset;
        n_slices = n_slices,
        tablica_k = tablica_k,
        feature_cols = [2, 3, 4],
        methods = methods
    )

    # 5. Generate individual plots for each regularization method
    output_dir = joinpath(@__DIR__, "..", "..", "plots")
    mkpath(output_dir)

    save(joinpath(output_dir, "hjsw_lpca_regularizations_comparison.png"), fig_comparison)

    for m in methods
        fig_m = plot_local_pca(
            dataset;
            n_slices = n_slices,
            tablica_k = tablica_k,
            feature_cols = [2, 3, 4],
            normalize = m,
            title = "LPCA z regularyzacją: $(m)"
        )
        save(joinpath(output_dir, "hjsw_lpca_method_$(m).png"), fig_m)
    end

    # 6. Generate Phase Space plot with LPCA dimensions colored per point
    println("Generating 2D Phase Space plot with LPCA dimensions colored per point...")
    fig_ps_2d = plot_phase_space_lpca_dims(
        dataset;
        liczba_sąsiadów = 15,
        wybrane_czasy = 0.25,
        metoda_regularyzacji = :max,
        zmienna_x = :T,
        zmienna_y = :A,
        tytuł_wykresu = "Wymiar LPCA w przestrzeni fazowej (T, A) dla τ = 0.25 fm/c"
    )
    save(joinpath(output_dir, "hjsw_phase_space_lpca_dims_2d.png"), fig_ps_2d)

    println("Generating 3D Phase Space plot with LPCA dimensions colored per point...")
    fig_ps_3d = plot_phase_space_lpca_dims(
        dataset;
        liczba_sąsiadów = 15,
        wybrane_czasy = [0.25, 0.5, 1.0],
        metoda_regularyzacji = :max,
        zmienna_x = :T,
        zmienna_y = :A,
        zmienna_z = :B,
        tytuł_wykresu = "Wymiar LPCA w przestrzeni fazowej 3D (T, A, B)"
    )
    save(joinpath(output_dir, "hjsw_phase_space_lpca_dims_3d.png"), fig_ps_3d)

    println("All regularization & phase space LPCA plots saved to: ", output_dir)
    return (fig_comparison = fig_comparison, fig_ps_2d = fig_ps_2d, fig_ps_3d = fig_ps_3d)
end

if abspath(PROGRAM_FILE) == @__FILE__
    run_hjsw_lpca_regularization_example()
end
