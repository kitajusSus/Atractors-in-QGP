using AttractorsQGP
using GLMakie
using LaTeXStrings

# Set working directory to the project root or resolve paths relative to this script
script_dir = @__DIR__
dataset_path = joinpath(script_dir, "../../datasets/ncbj_xin_2020_5000.h5")

if !isfile(dataset_path)
    println("Dataset not found at: ", dataset_path)
    println("Looking for fallback dataset...")
    dataset_path = joinpath(script_dir, "../../datasets/dataset_small_testowy_1000_mis.h5")
end

if !isfile(dataset_path)
    error("No datasets found to plot!")
end

println("Loading dataset: ", dataset_path)
dataset = load_dataset(dataset_path)

println("Generating plot...")
# We compute with 20 slices to balance performance and visual detail
fig = plot_lpca_entropy(dataset; zakres_K = [10, 20, 40, 80], n_slices = 20)

output_path = joinpath(script_dir, "plot_lpca_entropy.png")
println("Saving plot to: ", output_path)
save(output_path, fig)
println("Plot saved successfully!")
