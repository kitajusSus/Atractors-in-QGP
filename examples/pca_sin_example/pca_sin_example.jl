using GLMakie
using Random
using LinearAlgebra
using Statistics
using AttractorsQGP

function generate_noisy_sine(N; x_min = 0.0, x_max = 2 * pi, noise_std = 0.3, seed = 5)
    Random.seed!(seed)
    x = range(x_min, x_max, length = N) |> collect
    y_true = sin.(x)
    y_noisy = y_true .+ noise_std .* randn(N) * 0.1

    return hcat(x, y_noisy), x, y_true
end

function run_pca(X)
    X_norm, mn, mx = normalize_minmax(X)

    cov_matrix = cov(X_norm)

    eigen_decomp = eigen(cov_matrix)

    idx = sortperm(eigen_decomp.values, rev = true)
    values = eigen_decomp.values[idx]
    vectors = eigen_decomp.vectors[:, idx] # kolumny to wektory własne (PC1, PC2)

    mean_X_norm = mean(X_norm, dims = 1)

    return X_norm, mean_X_norm, values, vectors
end

function main()
    N_points = 250
    noise_level = 0.25

    println("Generowanie danych sinusa...")
    X, x_true, y_true = generate_noisy_sine(N_points; noise_std = noise_level)

    println("Obliczanie PCA...")
    X_norm, mean_X_norm, eigenvalues, eigenvectors = run_pca(X)

    println("\nWyniki PCA:")
    println("Średnia znormalizowanych danych (środek): ", mean_X_norm)
    println("Wartości własne (wariancja składowych): ", eigenvalues)
    println("Pierwsza główna składowa (PC1): ", eigenvectors[:, 1])
    println("Druga główna składowa (PC2): ", eigenvectors[:, 2])

    println("\nAplikowanie motywu set_publication_theme()...")
    set_publication_theme(cmap = :davos)

    fig = Figure(size = (1400, 650))

    ax1 = Axis(
        fig[1, 1],
        title = "Oryginalne dane zaszumione",
        xlabel = "x",
        ylabel = "y"
    )

    lines!(ax1, x_true, y_true, color = :crimson, linewidth = 3, label = "Czysty sinus")
    scatter!(ax1, X[:, 1], X[:, 2], color = RGBAf(0.1176, 0.5647, 1.0, 0.5), markersize = 10, label = "Zaszumione dane")
    axislegend(ax1, position = :rt)

    ax2 = Axis(
        fig[1, 2],
        title = "Normalizacja min-max i składowe główne (PCA)",
        xlabel = "x (znormalizowane)",
        ylabel = "y (znormalizowane)"
    )

    scatter!(ax2, X_norm[:, 1], X_norm[:, 2], color = RGBAf(0.1176, 0.5647, 1.0, 0.4), markersize = 10, label = "Dane znormalizowane")

    scale_factor = 2.0
    len1 = scale_factor * sqrt(eigenvalues[1])
    len2 = scale_factor * sqrt(eigenvalues[2])

    pc1_dir = eigenvectors[:, 1]
    pc2_dir = eigenvectors[:, 2]

    arrows!(
        ax2, [mean_X_norm[1]], [mean_X_norm[2]], [len1 * pc1_dir[1]], [len1 * pc1_dir[2]],
        color = :crimson, linewidth = 4, arrowsize = 22
    )
    arrows!(
        ax2, [mean_X_norm[1]], [mean_X_norm[2]], [len2 * pc2_dir[1]], [len2 * pc2_dir[2]],
        color = :forestgreen, linewidth = 4, arrowsize = 22
    )

    elem_data = MarkerElement(color = RGBAf(0.1176, 0.5647, 1.0, 0.4), marker = :circle, markersize = 12, strokecolor = :black, strokewidth = 0.5)
    elem_pc1 = LineElement(color = :crimson, linewidth = 4)
    elem_pc2 = LineElement(color = :forestgreen, linewidth = 4)

    axislegend(ax2, [elem_data, elem_pc1, elem_pc2], ["Dane znormalizowane", "PC1 ", "PC2 "], position = :rt)

    combined_path = joinpath(@__DIR__, "pca_sinus_combined.pdf")
    println("Zapisywanie połączonego wykresu do: ", combined_path)
    save(combined_path, fig)

    fig_data = Figure(size = (800, 600))
    ax_data = Axis(fig_data[1, 1], title = "Oryginalne dane zaszumione", xlabel = "x", ylabel = "y")
    lines!(ax_data, x_true, y_true, color = :crimson, linewidth = 3, label = "Czysty sinus")
    scatter!(ax_data, X[:, 1], X[:, 2], color = RGBAf(0.1176, 0.5647, 1.0, 0.5), markersize = 10, label = "Zaszumione dane")
    axislegend(ax_data, position = :rt)
    data_path = joinpath(@__DIR__, "noisy_sine_data.pdf")
    println("Zapisywanie wykresu danych do: ", data_path)
    save(data_path, fig_data)

    fig_pca = Figure(size = (800, 600))
    ax_pca = Axis(fig_pca[1, 1], title = "Analiza PCA na zaszumionym sinusie", xlabel = "x (znormalizowane)", ylabel = "y (znormalizowane)")
    scatter!(ax_pca, X_norm[:, 1], X_norm[:, 2], color = RGBAf(0.1176, 0.5647, 1.0, 0.4), markersize = 10, label = "Dane znormalizowane")
    arrows!(
        ax_pca, [mean_X_norm[1]], [mean_X_norm[2]], [len1 * pc1_dir[1]], [len1 * pc1_dir[2]],
        color = :crimson, linewidth = 4, arrowsize = 22
    )
    arrows!(
        ax_pca, [mean_X_norm[1]], [mean_X_norm[2]], [len2 * pc2_dir[1]], [len2 * pc2_dir[2]],
        color = :forestgreen, linewidth = 4, arrowsize = 22
    )
    axislegend(ax_pca, [elem_data, elem_pc1, elem_pc2], ["Dane znormalizowane", "PC1 ", "PC2 "], position = :rt)
    pca_path = joinpath(@__DIR__, "pca_arrows.pdf")
    println("Zapisywanie wykresu PCA do: ", pca_path)
    save(pca_path, fig_pca)

    return println("Gotowe!")
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
