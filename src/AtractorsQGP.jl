module AtractorsQGP


abstract type AbstractHydroModel end
"""
Hydrodynamic transport parameters.
"""
struct HydroParams{T <: Real}
    eta_over_s::T
    tau_pi::T
    lambda1::T
end

include("models/brsss.jl")
include("models/mis.jl")
include("models/hjsw.jl")

include("constants/units.jl")

include("equations/bjorken.jl")


# include("examples/lagranz.jl")
# export lle_ncbj, fun_y
# # funkcje do obliczen symbolicznych
# export lagrang, basic_lag
include("solver/hydro_solver.jl")

include("simulation/initial_conditions.jl")
include("simulation/trajectories.jl")

include("analysis/lle.jl")
export run_lle_per_time, lle, run_lle_for_selected_taus, lle_spectrum, lle_spectrum_over_k, spectrum_statistics, scan_lle_spectrum
include("analysis/lpca.jl")
export compute_lpca, normalizuj_max, dynamic_lpca_analysis, compute_lpca_entropy, compute_stable_lpca_collapse, compute_lpca_principal_angles
include("analysis/pca.jl")
include("analysis/dimension.jl")
include("analysis/plots.jl")
export animate_pca_evolution, plot_pca_bar_variance, plot_lle_grid, plot_lle_grid_2d, plot_lle_grid_3d, plot_lle_embedding
export plot_lle_results_for_taus, plot_lle_spectrum_analysis, plot_lle_spectrum_scan_analysis
export plot_lle_spectrum_statistics_grid, plot_lle_σ_spectrum
export plot_pinn_deff_evolution
export plot_phase_space_evolution, plot_phase_space_evolution_3d
export plot_local_pca
include("analysis/fit_polynomials.jl")
export compute_polynomial_lle

# ── PINNs ──────────────────────────────────────────────────────────────────
# include("pinns/network.jl")
# include("pinns/losses.jl")
# include("pinns/training.jl")
# include("pinns/pinn_solver.jl")
# include("pinns/jacobian_analysis.jl")
# include("pinns/data_jacobian_analysis.jl")
# export PINNConfig, build_pinn_network, normalize_pinn_input, denormalize_pinn_output
# export pinn_predict
# export train_pinn, PINNResult
# export predict_trajectories_pinn, build_pinn_dataset
# export compare_pinn_ode, pinn_attractor_analysis
# # Jacobian-based dimensionality reduction
# export pinn_jacobian, pinn_jacobian_full
# export d_eff_from_singular_values, pinn_deff_at
# export pinn_deff_scan, pinn_jacobian_scan
# export pinn_hydrodynamisation_time
# export sample_ic_ensemble, fixed_transport_ic_ensemble
# export pinn_dimensionality_workflow
# # Data-driven Neural Jacobians
# export normalize_generic, build_generic_network, compute_data_jacobian
include("io/data_io.jl")
include("io/attractor.jl")
export build_attractor_interpoland, get_attractor_line, attractor_data
# idl twonn
include("analysis/intrinsic_dimension.jl")
export estimate_lid, estimate_twonn, estimate_dimension, scan_intrinsic_dimensions


export HydroParams, AbstractHydroModel, BRSSSModel, MISModel, HJSWModel
export HBARC_MEV_FM, MEV_PER_FM, FM_PER_MEV, to_temperature_unit, temperature_to_fm
export solve_hydro, generate_initial_conditions, generate_trajectories
export build_dataset, run_pca, run_pca_kernel, run_pca_per_time, estimate_dimension, scan_dimension_from_data
export explained_variance_ratio_from_svd,
    normalize_minmax,
    run_pca,
    run_pca_kernel,
    get_tau_slice,
    run_pca_for_tau,
    run_pca_per_time,
    run_evolution_pca_workflow
export lle, plot_lle_dim, plot_simulation_lle, plot_lle_spectrum_statistics, plot_lle_spectrum_statistics_grid
export save_dataset_csv, load_dataset_csv, save_dataset_h5, load_dataset_h5
export save_dataset_jls, load_dataset_jls, save_dataset, load_dataset
export set_publication_theme, resolve_def, get_data, get_limits, get_attractor_line_for_frame
export plot_phase_space_grid, plot_thermodynamics_evolution, plot_pca_summary, plot_pca_evr_over_time, plot_twonn, plot_lid_dimension
export plot_map_lpca, plot_attractor_Aw_T, plot_attractor, plot_lpca_entropy, plot_stable_lpca_collapse, plot_lpca_principal_angles
export run_main


# include("../examples/ncbj_lle.jl")
# export ncbj1_macierz_wszystkich_punktow, ncbj2_sasiedzi, ncbj3_calculate_wagi_dla_x_i
# export ncbj3_svd_wagi_dla_x_i, ncbj4_lle_basic, ncbj4_lle_svd, ncbj5_nowy_manifold
# include("examples/lle_examples/plots_lle.jl")
# include("examples/lle_examples/swissroll.jl")
# include("examples/lle_examples/scurve.jl")
# include("examples/lle_examples/helix.jl")
#
# export ex_swissroll, swissroll_dane
# export ex_scurve, scurve_dane
# export ex_helix, helix_dane
#

"""
Run full simulation generating important data
"""
function run_main(
        model::AbstractHydroModel;
        n_points::Integer = 5000,
        tspan::Tuple{<:Real, <:Real} = (0.2, 1.2),
        T_range::Tuple{<:Real, <:Real} = (200.0, 1400.0),
        A_range::Tuple{<:Real, <:Real} = (-1, 7),
        saveat::Union{Real, AbstractVector{<:Real}, Nothing} = 0.01,
        seed::Integer = 5,
    )
    ics = generate_initial_conditions(n_points; T_range = T_range, A_range = A_range, seed = seed)
    solutions = generate_trajectories(model, ics, tspan; saveat = saveat)
    dataset = build_dataset(solutions)
    return (solutions = solutions, dataset = dataset)
end


function run_main(
        model::HJSWModel;
        n_points::Integer = 5000,
        tspan::Tuple{<:Real, <:Real} = (0.2, 1.2),
        T_range::Tuple{<:Real, <:Real} = (200.0, 1400.0),
        A_MIS_range::Union{Nothing, Tuple{<:Real, <:Real}} = nothing,
        A_QNM_range::Tuple{<:Real, <:Real} = (-2.0, 2.0),
        A_range::Union{Nothing, Tuple{<:Real, <:Real}} = nothing,
        saveat::Union{Real, AbstractVector{<:Real}, Nothing} = 0.01,
        seed::Integer = 5,
    )
    actual_A_MIS_range = if !isnothing(A_MIS_range)
        A_MIS_range
    elseif !isnothing(A_range)
        A_range
    else
        (-1.0, 10.0)
    end
    ics = generate_initial_conditions(
        model,
        n_points;
        T_range = T_range,
        A_MIS_range = actual_A_MIS_range,
        A_QNM_range = A_QNM_range,
        seed = seed,
    )
    solutions = generate_trajectories(model, ics, tspan; saveat = saveat)
    dataset = build_dataset(solutions)
    return (solutions = solutions, dataset = dataset)
end


end
