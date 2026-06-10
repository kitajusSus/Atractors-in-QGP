using Test
# include("../src/AtractorsQGP.jl")
using AtractorsQGP
using StaticArrays
using Random
using Lux
using ComponentArrays
using GLMakie

@testset "AtractorsQGP refactor API" begin
    model = BRSSSModel()
    u0 = [800.0, 0.5]
    tspan = (0.22, 0.5)

    sol = solve_hydro(model, u0, tspan; saveat = 0.02)
    @test !isempty(sol.t)
    τ = 0.3
    T, A = 5.0, 0.4
    p = model.params
    du = AtractorsQGP.rhs([T, A], model, τ)
    expected_dT = (T / τ) * (-1 / 3 + A / 18)
    expected_term_T = τ * T * (A + (p.lambda1 / (12 * p.eta_over_s)) * A^2)
    expected_term_A2 = (2 / 9) * p.tau_pi * A^2
    expected_dA = (1 / (p.tau_pi * τ)) * (8 * p.eta_over_s - expected_term_T - expected_term_A2)
    @test du[1] ≈ expected_dT atol = 1.0e-12
    @test du[2] ≈ expected_dA atol = 1.0e-12


    ics_default_a = generate_initial_conditions(5; T_range = (700.0, 900.0), A_range = (-1.0, 1.0))
    ics_default_b = generate_initial_conditions(5; T_range = (700.0, 900.0), A_range = (-1.0, 1.0))
    @test ics_default_a == ics_default_b

    ics = ics_default_a
    @test eltype(ics) == SVector{2, Float64}
    @test all(ic -> 700.0 * FM_PER_MEV <= ic[1] <= 900.0 * FM_PER_MEV, ics)

    ics_fm = generate_initial_conditions(5; T_range = (3.0, 4.0), temperature_unit = :fm, rng = MersenneTwister(123))
    @test all(ic -> 3.0 <= ic[1] <= 4.0, ics_fm)

    sols = generate_trajectories(model, ics, tspan; saveat = 0.02, parallel = :serial)
    @test length(sols) == 5
    @test eltype(sols) != Any

    dataset = build_dataset(sols)
    @test size(dataset, 2) == 3

    dataset_mev = build_dataset(sols; temperature_unit = :MeV)
    @test size(dataset_mev) == size(dataset)
    @test dataset_mev[1, 2] ≈ dataset[1, 2] * MEV_PER_FM
    @test_throws ArgumentError build_dataset(sols; temperature_unit = :Kelvin)

    split_input = [
        0.2 1000.0 1.0
        0.3  900.0 0.5
        0.2  850.0 0.2
        0.3  800.0 0.1
    ]
    split_ranges = AtractorsQGP._split_trajectories(split_input)
    @test split_ranges == [1:2, 3:4]

    bad_solutions = [
        (t = [0.2, 0.21], u = [SVector(1000.0, 0.5), SVector(NaN, 0.4)]),
        (t = [0.2], u = [SVector(950.0, -0.2)]),
    ]
    cleaned = build_dataset(bad_solutions)
    @test size(cleaned, 1) == 2
    @test all(isfinite, cleaned)

    # lle = run_LLE(model, u0, tspan; saveat=0.02)
    # @test isfinite(lle)

    dim = estimate_dimension(dataset[:, 2:3])
    @test dim > 0
    flow = run_evolution_pca_workflow(
        model;
        n_points = 6,
        tspan = (0.22, 0.3),
        saveat = 0.02,
        parallel = :serial,
        feature_cols = [2, 3],
        rng = MersenneTwister(7),
    )
    @test length(flow.initial_conditions) == 6
    @test length(flow.solutions) == 6
    @test size(flow.dataset, 2) == 3
    @test size(flow.pca_over_time.explained_variance_ratio, 2) == 2
    flow_default_a = run_evolution_pca_workflow(
        model;
        n_points = 4,
        tspan = (0.22, 0.28),
        saveat = 0.02,
        parallel = :serial,
    )
    flow_default_b = run_evolution_pca_workflow(
        model;
        n_points = 4,
        tspan = (0.22, 0.28),
        saveat = 0.02,
        parallel = :serial,
    )
    @test flow_default_a.dataset == flow_default_b.dataset


    # PCA per-time should honor selected feature columns
    pca_dataset = [
        0.2  1.0   10.0
        0.2  2.0   20.0
        0.2  3.0   30.0
        0.3  10.0   1.0
        0.3  20.0   2.0
        0.3  30.0   3.0
    ]
    tau_only = sort(Float64.(pca_dataset[:, 1]))
    idx_tau = get_tau_slice(tau_only, 0.2)
    @test length(idx_tau) == 3

    _, x_default = get_tau_slice(pca_dataset, 0.2)
    _, x_selected = get_tau_slice(pca_dataset, 0.2; feature_cols = [3])
    @test size(x_default, 2) == 2
    @test size(x_selected, 2) == 1
    @test x_selected[:, 1] == [10.0, 20.0, 30.0]

    pca_time = run_pca_per_time(pca_dataset; n_components = 1, feature_cols = [2])
    @test size(pca_time.explained_variance_ratio) == (2, 1)
    @test all(isfinite, pca_time.explained_variance_ratio)
    one_tau = run_pca_for_tau(pca_dataset, 0.2; n_components = 1, feature_cols = [2, 3])
    @test one_tau.n_points == 3
    @test length(one_tau.pca_result.explained_variance_ratio) == 1


    shifted = [1.0 1.0; 2.0 2.0; 3.0 3.0]
    shifted_plus = shifted .+ 100.0
    evr_a = run_pca(shifted; n_components = 2).explained_variance_ratio
    evr_b = run_pca(shifted_plus; n_components = 2).explained_variance_ratio
    @test evr_a ≈ evr_b atol = 1.0e-12


    # PCA per-time should honor selected feature columns
    pca_dataset = [
        0.2  1.0   10.0
        0.2  2.0   20.0
        0.2  3.0   30.0
        0.3  10.0   1.0
        0.3  20.0   2.0
        0.3  30.0   3.0
    ]
    tau_only = sort(Float64.(pca_dataset[:, 1]))
    idx_tau = get_tau_slice(tau_only, 0.2)
    @test length(idx_tau) == 3

    _, x_default = get_tau_slice(pca_dataset, 0.2)
    _, x_selected = get_tau_slice(pca_dataset, 0.2; feature_cols = [3])
    @test size(x_default, 2) == 2
    @test size(x_selected, 2) == 1
    @test x_selected[:, 1] == [10.0, 20.0, 30.0]

    pca_time = run_pca_per_time(pca_dataset; n_components = 1, feature_cols = [2])
    @test size(pca_time.explained_variance_ratio) == (2, 1)
    @test all(isfinite, pca_time.explained_variance_ratio)
    one_tau = run_pca_for_tau(pca_dataset, 0.2; n_components = 1, feature_cols = [2, 3])
    @test one_tau.n_points == 3
    @test length(one_tau.pca_result.explained_variance_ratio) == 1


    shifted = [1.0 1.0; 2.0 2.0; 3.0 3.0]
    shifted_plus = shifted .+ 100.0
    evr_a = run_pca(shifted; n_components = 2).explained_variance_ratio
    evr_b = run_pca(shifted_plus; n_components = 2).explained_variance_ratio
    @test evr_a ≈ evr_b atol = 1.0e-12


    mktempdir() do tmp
        csv_path = joinpath(tmp, "dataset.csv")
        h5_path = joinpath(tmp, "dataset.h5")
        jls_path = joinpath(tmp, "dataset.jls")

        save_dataset(csv_path, dataset)
        save_dataset(h5_path, dataset)
        save_dataset(jls_path, dataset)

        @test size(load_dataset(csv_path)) == size(dataset)
        @test size(load_dataset(h5_path)) == size(dataset)
        @test size(load_dataset(jls_path)) == size(dataset)

        csv_lower = joinpath(tmp, "dataset_lower.csv")
        write(csv_lower, "tau,t,a
0.22,1.0,2.0
0.23,1.1,2.1
")
        lower_loaded = load_dataset(csv_lower)
        @test size(lower_loaded) == (2, 3)
    end

    @testset "get_limits and plot_phase_space_grid limits" begin
        # sample dataset
        data = [
            0.2  1.0  -1.0
            0.2  2.0   1.0
            0.3  1.5  -0.5
            0.3  2.5   0.5
        ]
        
        # Test single axis limits
        xlims = get_limits(data, :T; times=(0.2, 0.3))
        @test xlims[1] < 1.0
        @test xlims[2] > 2.5
        
        ylims = get_limits(data, :A; times=(0.2, 0.3))
        @test ylims[1] < -1.0
        @test ylims[2] > 1.0

        # Test both axes limits
        lims = get_limits(data, :T, :A; times=(0.2, 0.3))
        @test lims[1] == xlims[1]
        @test lims[2] == xlims[2]
        @test lims[3] == ylims[1]
        @test lims[4] == ylims[2]
    end

    @testset "scan_dimension_from_data" begin
        mock_data = [
            0.2  1.0  -1.0  0.5
            0.2  2.0   1.0  0.6
            0.2  3.0   0.0  0.4
            0.2  1.5  -0.5  0.5
            0.3  1.5  -0.5  1.0
            0.3  2.5   0.5  2.0
            0.3  2.0   0.0  1.5
            0.3  1.8  -0.2  1.2
        ]
        scan_all = scan_dimension_from_data(mock_data)
        @test scan_all.taus == [0.2, 0.3]
        @test length(scan_all.dimensions) == 2
        @test all(d -> 1.0 <= d <= 3.0, scan_all.dimensions)

        scan_subset = scan_dimension_from_data(mock_data; feature_cols=[2, 3])
        @test scan_subset.taus == [0.2, 0.3]
        @test length(scan_subset.dimensions) == 2
        @test all(d -> 1.0 <= d <= 2.0, scan_subset.dimensions)
    end

    @testset "Data-driven Neural Jacobians" begin
        mock_data = [
            0.1 10.0 100.0
            0.5 50.0 500.0
        ]
        norm_data, mn, mx = normalize_generic(mock_data)
        @test size(norm_data) == (2, 3)
        @test all(v -> -1.0 <= v <= 1.0, norm_data)
        @test mn ≈ [0.1, 10.0, 100.0]
        @test mx ≈ [0.5, 50.0, 500.0]

        nn = build_generic_network(3, 2; hidden_layers=2, hidden_size=10)
        rng = Random.MersenneTwister(42)
        ps_raw, st = Lux.setup(rng, nn)
        ps = ComponentArray{Float64}(ps_raw)

        mn_in = [0.1, 10.0, 100.0]
        mx_in = [0.5, 50.0, 500.0]
        mn_out = [1.0, 10.0]
        mx_out = [2.0, 20.0]

        τ = 0.3
        x0 = [30.0, 300.0]

        J = compute_data_jacobian(nn, ps, st, τ, x0, mn_in, mx_in, mn_out, mx_out)
        @test size(J) == (2, 2)
        @test all(isfinite, J)

        # 3. Test plot_pinn_deff_evolution
        wyniki_mock = [
            0.2 1.5
            0.3 1.2
            0.4 1.05
        ]
        fig = plot_pinn_deff_evolution(wyniki_mock)
        @test fig isa GLMakie.Figure

        # 4. Test plot_lle_spectrum_statistics_grid
        mock_model = BRSSSModel()
        mock_ics = generate_initial_conditions(15; T_range=(700.0, 900.0), A_range=(-1.0, 1.0))
        mock_sols = generate_trajectories(mock_model, mock_ics, (0.22, 0.3); saveat=0.04, parallel=:serial)
        mock_dataset = build_dataset(mock_sols)
        fig_grid = plot_lle_spectrum_statistics_grid(mock_dataset, [0.22, 0.26]; k_values = 2:2:4, ile_λ = 2)
        @test fig_grid isa GLMakie.Figure
    end

    @testset "HJSW Model ewolucja i PCA" begin
        model = HJSWModel()
        @test model.omega_R ≈ 9.8
        @test model.omega_I ≈ 8.629

        ics = generate_initial_conditions(model, 5; T_range=(800.0, 1000.0), A_MIS_range=(-1.0, 1.0), A_QNM_range=(-0.5, 0.5))
        @test eltype(ics) == SVector{4, Float64}
        @test length(ics) == 5

        sols = generate_trajectories(model, ics, (0.22, 0.3); saveat=0.02, parallel=:serial)
        @test length(sols) == 5

        dataset = build_dataset(sols)
        @test size(dataset, 2) == 5 # [tau, T, A_MIS, A_QNM, dA_QNM]
        @test all(isfinite, dataset)

        # check that PCA and dimension scanner work on this 5-column dataset:
        wynik_dim = scan_dimension_from_data(dataset)
        @test length(wynik_dim.taus) > 0
        @test all(d -> 1.0 <= d <= 4.0, wynik_dim.dimensions)

        wynik_pca = run_pca_per_time(dataset; n_components=2)
        @test length(wynik_pca.taus) > 0
        @test size(wynik_pca.explained_variance_ratio, 2) == 2

        # test run_main with HJSWModel
        hjsw_res = run_main(model, n_points = 5, tspan = (0.2, 0.3), T_range = (200.0, 1400.0), A_range = (-1.0, 8.0), seed = 5)
        @test length(hjsw_res.solutions) == 5
        @test size(hjsw_res.dataset, 2) == 5
    end
end
