"""
    plot_lpca_entropy(
        dataset::AbstractMatrix{<:Real};
        zakres_K::AbstractVector{<:Integer} = [10, 20, 40, 80],
        n_slices::Union{Nothing, Int} = nothing,
        feature_cols::AbstractVector{<:Integer} = collect(2:size(dataset, 2))
    )

Plots the mean local PCA entropy (z odchyleniem standardowym) as a function of time tau
for different values of K.
"""
function plot_lpca_entropy(
        dataset::AbstractMatrix{<:Real};
        zakres_K::AbstractVector{<:Integer} = [10, 20, 40, 80],
        n_slices::Union{Nothing, Int} = nothing,
        feature_cols::AbstractVector{<:Integer} = collect(2:size(dataset, 2))
    )

    set_publication_theme()

    fig = Figure(size = (950, 600))
    ax = Axis(
        fig[1, 1],
        title = L"\text{Średnia lokalna entropia spektralna PCA } \langle S_{\mathrm{PCA}} \rangle \text{ jako funkcja } \tau",
        xlabel = L"\tau\,[\mathrm{fm}/c]",
        ylabel = L"\langle S_{\mathrm{PCA}} \rangle",
        xticks = 0:0.25:maximum(dataset[:, 1]),
        limits = (nothing, nothing, -0.05, 1.05),
        xautolimitmargin = (0.02, 0.02),
        yautolimitmargin = (0.05, 0.05)
    )
    palette = Makie.theme(:Palette).color[]
    # palette = [:crimson, :dodgerblue, :forestgreen, :darkorange, :purple]

    for (i, k) in enumerate(zakres_K)
        tau_vals, mean_entropy, std_entropy = compute_lpca_entropy(
            dataset, k; n_slices = n_slices, feature_cols = feature_cols
        )

        c = palette[mod1(i, length(palette))]
        band!(
            ax,
            tau_vals,
            max.(mean_entropy .- std_entropy, 0.0),
            min.(mean_entropy .+ std_entropy, 1.0);
            color = (c, 0.1)
        )

        lines!(
            ax,
            tau_vals,
            mean_entropy;
            linewidth = 3.0,
            color = palette[mod1(i * 2, length(palette))],
            label = L"K = %$(k)"
        )
    end

    axislegend(ax, position = :rt)

    return fig
end

"""
    plot_stable_lpca_collapse(
        dataset::AbstractMatrix{<:Real};
        zakres_K::AbstractVector{<:Integer} = [10, 20, 30, 40],
        Scrit::Real = 0.2,
        delta_k::Real = 0.05,
        n_slices::Union{Nothing, Int} = 30,
        feature_cols::AbstractVector{<:Integer} = collect(2:size(dataset, 2)),
        norm_type::Symbol = :embedding
    )

Computes and plots the stable LPCA collapse analysis.
Creates a 2-panel figure showing:
1. The mean spectral entropy <S_PCA> over scales k with a band for the spread.
2. The scale spread delta_k S_PCA.
Vertical lines mark the stable collapse time tau_LPCA.
"""
function plot_stable_lpca_collapse(
        dataset::AbstractMatrix{<:Real};
        zakres_K::AbstractVector{<:Integer} = [10, 20, 30, 40],
        Scrit::Real = 0.2,
        delta_k::Real = 0.05,
        n_slices::Union{Nothing, Int} = 30,
        feature_cols::AbstractVector{<:Integer} = collect(2:size(dataset, 2)),
        norm_type::Symbol = :embedding
    )
    set_publication_theme()

    # Compute stable collapse data
    res = compute_stable_lpca_collapse(
        dataset;
        zakres_K = zakres_K,
        Scrit = Scrit,
        delta_k = delta_k,
        n_slices = n_slices,
        feature_cols = feature_cols,
        norm_type = norm_type
    )

    fig = Figure(size = (1000, 750))

    # Panel 1: Mean Entropy
    ax1 = Axis(
        fig[1, 1],
        title = L"\text{Stabilny kolaps wymiaru LPCA: } \langle S_{\mathrm{PCA}} \rangle_k",
        ylabel = L"\langle S_{\mathrm{PCA}} \rangle_k",
        xticks = 0:0.25:maximum(res.taus),
        limits = (nothing, nothing, -0.05, 1.05)
    )

    # Panel 2: Spread over scales
    ax2 = Axis(
        fig[2, 1],
        xlabel = L"\tau\,[\mathrm{fm}/c]",
        ylabel = L"\Delta_k S_{\mathrm{PCA}}",
        xticks = 0:0.25:maximum(res.taus),
        limits = (nothing, nothing, -0.01, maximum(res.S_PCA_std) * 1.3 + 0.05)
    )

    # Plot on panel 1
    # Plot individual curves for each k
    palette = [:crimson, :dodgerblue, :forestgreen, :darkorange, :purple]
    for (i, k_val) in enumerate(zakres_K)
        c = palette[mod1(i, length(palette))]
        lines!(
            ax1,
            res.taus,
            res.S_PCA_matrix[:, i];
            color = (c, 0.45),
            linewidth = 1.5,
            linestyle = :dot,
            label = L"K = %$(k_val)"
        )
    end

    # Plot the mean and band
    band!(
        ax1,
        res.taus,
        max.(res.S_PCA_mean .- res.S_PCA_std, 0.0),
        min.(res.S_PCA_mean .+ res.S_PCA_std, 1.0);
        color = (:black, 0.15)
    )
    lines!(
        ax1,
        res.taus,
        res.S_PCA_mean;
        color = :black,
        linewidth = 3.5,
        label = L"\text{Średnia } \langle S_{\mathrm{PCA}} \rangle_k"
    )

    hlines!(ax1, [Scrit]; color = :red, linestyle = :dash, linewidth = 2.0, label = L"S_{\mathrm{crit}}")

    lines!(
        ax2,
        res.taus,
        res.S_PCA_std;
        color = :blue,
        linewidth = 3.0
    )
    # Threshold delta_k
    hlines!(ax2, [delta_k]; color = :red, linestyle = :dash, linewidth = 2.0, label = L"\delta_k")

    # If stable collapse is found, draw vertical lines
    if !isnan(res.tau_LPCA_stable)
        vlines!(ax1, [res.tau_LPCA_stable]; color = :forestgreen, linestyle = :dashdot, linewidth = 3.0, label = L"\tau_{\mathrm{LPCA}} = %$(round(res.tau_LPCA_stable, digits=3))")
        vlines!(ax2, [res.tau_LPCA_stable]; color = :forestgreen, linestyle = :dashdot, linewidth = 3.0)
    end

    axislegend(ax1, position = :rt)
    axislegend(ax2, position = :rt)

    # Link x-axes for synchronized panning/zooming
    linkxaxes!(ax1, ax2)

    return fig, res
end


"""
    plot_lpca_principal_angles(
        dataset::AbstractMatrix{<:Real};
        k::Int = 20,
        subspace_dim::Int = 1,
        n_slices::Union{Nothing, Int} = 30,
        feature_cols::AbstractVector{<:Integer} = collect(2:size(dataset, 2))
    )

"""
function plot_lpca_principal_angles(
        dataset::AbstractMatrix{<:Real};
        k::Int = 20,
        subspace_dim::Int = 1,
        n_slices::Union{Nothing, Int} = 30,
        feature_cols::AbstractVector{<:Integer} = collect(2:size(dataset, 2))
    )
    set_publication_theme()

    res = compute_lpca_principal_angles(
        dataset;
        k = k,
        subspace_dim = subspace_dim,
        n_slices = n_slices,
        feature_cols = feature_cols
    )

    fig = Figure(size = (950, 600))
    ax = Axis(
        fig[1, 1],
        title = L"\text{Podobieństwo sąsiednich przestrzeni stycznych: kąty główne } \theta",
        xlabel = L"\tau\,[\mathrm{fm}/c]",
        ylabel = L"\text{Kąt główny } \theta\,[^\circ]",
        xticks = 0:0.25:maximum(res.taus),
        limits = (nothing, nothing, -5, 95)
    )

    # Plot band (mean ± std)
    band!(
        ax,
        res.taus,
        max.(res.mean_angles .- res.std_angles, 0.0),
        min.(res.mean_angles .+ res.std_angles, 90.0);
        color = (:dodgerblue, 0.15)
    )

    # Plot mean angle line
    lines!(
        ax,
        res.taus,
        res.mean_angles;
        color = :dodgerblue,
        linewidth = 3.5,
        label = L"\text{Średni kąt } \langle \theta \rangle"
    )

    # Draw a line for 90 degrees (maximum misalignment) and 0 degrees (perfect alignment)
    hlines!(ax, [0.0]; color = :gray, linestyle = :dot, linewidth = 1.0)
    hlines!(ax, [90.0]; color = :gray, linestyle = :dot, linewidth = 1.0)

    axislegend(ax, position = :rt)

    return fig, res
end
