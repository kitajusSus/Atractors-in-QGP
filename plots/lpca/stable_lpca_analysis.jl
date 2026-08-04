using AttractorsQGP
using GLMakie
using LaTeXStrings
using Statistics

# Set working directory relative to this script
script_dir = @__DIR__
dataset_path = joinpath(script_dir, "../../datasets/ncbj_xin_2020_5000.h5")

if !isfile(dataset_path)
    println("Dataset not found at: ", dataset_path)
    println("Looking for fallback dataset...")
    dataset_path = joinpath(script_dir, "../../datasets/dataset_small_testowy_1000_mis.h5")
end

if !isfile(dataset_path)
    error("No datasets found to run analysis!")
end

println("Loading dataset: ", dataset_path)
dataset = load_dataset(dataset_path)
println("Dataset size: ", size(dataset))

# -------------------------------------------------------------
# 1. Stable LPCA Collapse Analysis
# -------------------------------------------------------------
println("Running Stable LPCA Collapse analysis...")
zakres_K = [10, 20, 40, 80]
Scrit = 0.25
delta_k = 0.05
n_slices = 25

fig_collapse, res_collapse = plot_stable_lpca_collapse(
    dataset;
    zakres_K = zakres_K,
    Scrit = Scrit,
    delta_k = delta_k,
    n_slices = n_slices
)

println("Stable LPCA Collapse Results:")
println("  - Selected scales (K): ", zakres_K)
println("  - Entropy threshold (Scrit): ", Scrit)
println("  - Scale spread threshold (delta_k): ", delta_k)
for (i, k) in enumerate(zakres_K)
    println("  - Collapse time tau_LPCA(K=$(k)): ", res_collapse.tau_LPCA_k[i], " fm/c")
end
println("  - Stable collapse time tau_LPCA: ", res_collapse.tau_LPCA_stable, " fm/c")

collapse_plot_path = joinpath(script_dir, "stable_lpca_collapse.png")
println("Saving stable LPCA collapse plot to: ", collapse_plot_path)
save(collapse_plot_path, fig_collapse)

# -------------------------------------------------------------
# 2. Principal Angles between Local Tangent Spaces
# -------------------------------------------------------------
println("\nRunning Principal Angles analysis...")
k_neighbors = 20
subspace_dim = 1 # attractor is effectively 1D

fig_angles, res_angles = plot_lpca_principal_angles(
    dataset;
    k = k_neighbors,
    subspace_dim = subspace_dim,
    n_slices = n_slices
)

# Print some summary stats of the principal angles
println("Principal Angles Results (K=$(k_neighbors), d_subspace=$(subspace_dim)):")
println("  - First 5 time slices average angle:")
for i in 1:min(5, length(res_angles.taus))
    println("    tau = ", round(res_angles.taus[i], digits=3), " fm/c: mean angle = ", round(res_angles.mean_angles[i], digits=2), "° ± ", round(res_angles.std_angles[i], digits=2), "°")
end
println("  - Last 5 time slices average angle:")
n_t = length(res_angles.taus)
for i in max(1, n_t-4):n_t
    println("    tau = ", round(res_angles.taus[i], digits=3), " fm/c: mean angle = ", round(res_angles.mean_angles[i], digits=2), "° ± ", round(res_angles.std_angles[i], digits=2), "°")
end

angles_plot_path = joinpath(script_dir, "lpca_principal_angles.png")
println("Saving principal angles plot to: ", angles_plot_path)
save(angles_plot_path, fig_angles)

println("\nStable LPCA analysis completed successfully!")
