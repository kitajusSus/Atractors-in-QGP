using GLMakie 
using LaTeXStrings
using ColorSchemes
function set_publication_theme()
    scientific_palette = Makie.resample_cmap(:haline, 10) 

    set_theme!(
        Theme(
            font = "Libertinus Serif", 
            fontsize = 24, 
            figure_padding = 20,

            Axis = (
                backgroundcolor = :white,
                titlesize = 28,
                xlabelsize = 26,
                ylabelsize = 26,
                xticklabelsize = 20,
                yticklabelsize = 20,

                xgridstyle = :dash,
                ygridstyle = :dash,
                xgridcolor = RGBAf(0.85, 0.85, 0.85, 0.6),
                ygridcolor = RGBAf(0.85, 0.85, 0.85, 0.6),

                spinewidth = 1.5,
                bottomspinecolor = :black,
                leftspinecolor = :black,
                topspinecolor = :black,
                rightspinecolor = :black,
                topspinevisible = true,
                rightspinevisible = true,

               xtickalign = 1.0,
                ytickalign = 1.0,
                xticksize = 12,
                yticksize = 12,
                xtickwidth = 1.5,
                ytickwidth = 1.5,
                xtickcolor = :black,
                ytickcolor = :black,
                
                xminorticksvisible = true,
                yminorticksvisible = true,
                xminortickalign = 1.0,
                yminortickalign = 1.0,
                xminorticksize = 6,
                yminorticksize = 6,
                xminortickwidth = 1.0,
                yminortickwidth = 1.0,
            ),
            
            Legend = (
                framevisible = true,
                framewidth = 1.2,
                framecolor = :black,
                backgroundcolor = RGBAf(1.0, 1.0, 1.0, 0.90), 
                position = :rt,
                titlesize = 22,
                labelsize = 20,
                padding = (10.0, 10.0, 10.0, 10.0),
            ),

            Lines = (
                linewidth = 2.5, 
            ),
            Scatter = (
                markersize = 10,
                strokewidth = 0.5, 
                strokecolor = :black,
            ),

            Palette = ( 
                color = scientific_palette,
                patchcolor = scientific_palette, 
            ),
        )
    )
end
#
const PLOT_KEYS = Dict(
    :T => (L"T\,[\mathrm{fm}^{-1}]", x -> x[2]),
    :A => (L"\mathcal{A}", x -> x[3]),
    :tauT => (L"\tau T", x -> x[1] * x[2]),
    :tau2A => (L"\tau^2 \mathcal{A}", x -> x[1]^2 * x[3]),
)

function resolve_def(def)
    if def isa Symbol
        @assert haskey(PLOT_KEYS, def) "Unknown plot key: $def"
        return PLOT_KEYS[def]
    end

    if def isa Tuple && length(def) == 2
        return (def[1], def[2])
    end

    error("Axis definition must be Symbol or Tuple(label, function).")
end

function get_data(dataset::AbstractMatrix{<:Real}, t::Real, xdef, ydef)
    @assert size(dataset, 2) == 3 "Dataset must have columns [tau, T, A]."

    xlbl, xfn = resolve_def(xdef)
    ylbl, yfn = resolve_def(ydef)

    rows = findall(isapprox.(dataset[:, 1], t; atol = 1e-8))
    if isempty(rows)
        nearest = argmin(abs.(dataset[:, 1] .- t))
        rows = findall(isapprox.(dataset[:, 1], dataset[nearest, 1]; atol = 1e-8))
    end

    selected = dataset[rows, :]
    x = [xfn(selected[i, :]) for i = 1:size(selected, 1)]
    y = [yfn(selected[i, :]) for i = 1:size(selected, 1)]

    return (x = x, y = y, xlabel = xlbl, ylabel = ylbl)
end

function _split_trajectories(dataset::AbstractMatrix{<:Real})
    @assert size(dataset, 2) == 3 "Dataset must have columns [tau, T, A]."
    if size(dataset, 1) == 0
        return UnitRange{Int}[]
    end

    starts = Int[1]
    for i = 2:size(dataset, 1)
        if dataset[i, 1] <= dataset[i - 1, 1]
            push!(starts, i)
        end
    end

    ranges = UnitRange{Int}[]
    for k in eachindex(starts)
        s = starts[k]
        e = k < length(starts) ? starts[k + 1] - 1 : size(dataset, 1)
        push!(ranges, s:e)
    end
    return ranges
end

function plot_phase_space_grid(dataset::AbstractMatrix{<:Real}, times, xdef, ydef)
    set_publication_theme()
    
    palette = Makie.theme(:Palette).color[]
    n = length(times)
    ncols = min(3, n)
    nrows = ceil(Int, n / ncols)
    fig = Figure(size = (360 * ncols, 290 * nrows))

    for (i, t) in enumerate(times)
        row = (i - 1) ÷ ncols + 1
        col = (i - 1) % ncols + 1
        d = get_data(dataset, t, xdef, ydef)

        ax = Axis(
            fig[row, col],
            title = L"\tau = %$(round(t, digits=2))\,\mathrm{fm}/c",
            xlabel = d.xlabel,
            ylabel = d.ylabel,
        )
        
        scatter!(ax, d.x, d.y; markersize = 3.0, color = (palette[2], 0.80))
    end
    return fig
end

function plot_thermodynamics_evolution(dataset::AbstractMatrix{<:Real})
    set_publication_theme()

    palette = Makie.theme(:Palette).color[]
    trajs = _split_trajectories(dataset)
    fig = Figure(size = (950, 620))
    ax1 = Axis(
        fig[1, 1],
        title = L"\text{Ewolucja Temperatury } T\,[\mathrm{fm}^{-1}]\; \text{w czasie własnym } \tau",
        xlabel = L"\tau\,[\mathrm{fm}/c]",
        ylabel = L"T\,[\mathrm{fm}^{-1}]",
    )

    for tr in trajs
        lines!(
            ax1,
            dataset[tr, 1],
            dataset[tr, 2],
            color = (palette[1], 0.20),
            linewidth = 1.5,
        )
    end

    ax2 = Axis(
        fig[2, 1],
        title = L"\text{Ewolucja Anizotropii}\; \mathcal{A(τ)}\; \text{ w czasie własnym}",
        xlabel = L"\tau\,[\mathrm{fm}/c]",
        ylabel = L"\mathcal{A}",
    )

    for tr in trajs
        lines!(
            ax2,
            dataset[tr, 1],
            dataset[tr, 3],
            color = (palette[1], 0.20),
            linewidth = 1.5,
        )
    end

    hlines!(
        ax2,
        [0.0],
        color = palette[2], 
        linestyle = :dash,
        linewidth = 2.0,
        label = L"\mathcal{A}=0\;(\text{Anizotropia} = 0)",
    )
    axislegend(ax2, position = :rt)

    linkxaxes!(ax1, ax2)
    return fig
end

function plot_pca_evr_over_time(
    dataset::AbstractMatrix{<:Real};
    n_components::Int = 2,
    method::Symbol = :minmax,
    gamma::Float64 = 1.0,
    feature_cols::AbstractVector{<:Integer} = collect(2:size(dataset, 2)),
)
    set_publication_theme()

    palette = Makie.theme(:Palette).color[]
    result = run_pca_per_time(
        dataset;
        n_components = n_components,
        method = method,
        gamma = gamma,
        feature_cols = feature_cols,
    )
    taus = result.taus
    evr = result.explained_variance_ratio

    fig = Figure(size = (1600, 1000))
    ax = Axis(
        fig[1, 1],
        title = L"\text{Explained Variance (EVR) w funkcji czasu dla zmiennych } \mathcal{A} \text{ i } T/T_0",
        xlabel = L"\tau\,[\mathrm{fm}/c]",
        ylabel = L"\text{EVR}",
        limits = (0.20, maximum(taus), 0, 1.05),
    )

    hlines!(ax, [1.0], color = :gray45, linestyle = :dash, label = L"100\%")

    palette = Makie.theme(:Palette).color[]
    for comp = 1:n_components
        vals = evr[:, comp]
        mask = .!isnan.(vals)
        if any(mask)
            lines!(
                ax,
                taus[mask],
                vals[mask],
                linewidth = 3.0,
                color = palette[min(comp, length(palette))],
                label = L"\text{PC}%$(comp)",
            )
            if comp == 1
                band!(
                    ax,
                    taus[mask],
                    zeros(sum(mask)),
                    vals[mask],
                    color = (palette[1], 0.15),
                )
            end
        end
    end

    axislegend(ax, position = :rb)
    return fig
end

function plot_pca_summary(
    dataset::AbstractMatrix{<:Real};
    tau::Union{Nothing,Real} = nothing,
    tau_tol::Float64 = 1e-8,
    tau_mode::Symbol = :nearest,
    n_components::Int = 2,
    method::Symbol = :minmax,
    gamma::Float64 = 1.0,
)
    @assert size(dataset, 2) >= 3 "Dataset must contain at least [tau, T, A]."
    @assert tau_tol >= 0 "tau_tol must be >= 0."
    @assert tau_mode in (:strict, :nearest) "tau_mode must be :strict or :nearest."

    set_publication_theme()

    palette = Makie.theme(:Palette).color[]
    data = Matrix{Float64}(dataset)
    subtitle = L"\text{Wszystkie } \tau"

    if tau !== nothing
        τ = Float64(tau)
        τcol = data[:, 1]
        d = abs.(τcol .- τ)
        strict_mask = d .<= tau_tol

        if any(strict_mask)
            data = data[strict_mask, :]
            subtitle = L"\tau=%$(τ) \pm %$(tau_tol)"
        else
            if tau_mode === :strict
                error("No rows found for tau=$(τ) within tau_tol=$(tau_tol).")
            else
                i = argmin(d)
                τnearest = τcol[i]
                near_mask = τcol .== τnearest
                data = data[near_mask, :]
                subtitle = L"\text{najbliższe } \tau=%$(τnearest)"
            end
        end
    end

    features = data[:, 2:3]
    @assert size(features, 1) > 1 "Need at least two samples in selected tau slice."

    pca_result = if method === :minmax
        run_pca(features; n_components = n_components)
    elseif method === :kernel
        run_pca_kernel(features; n_components = n_components, gamma = gamma)
    else
        error("Unknown PCA method. Choose :minmax or :kernel.")
    end

    transformed = pca_result.transformed
    evr = pca_result.explained_variance_ratio
    n_show = min(size(transformed, 2), 2)

    fig = Figure(size = (980, 420))

    ax_proj = Axis(
        fig[1, 1],
        xlabel = L"\text{PC1}",
        ylabel = L"\text{PC2}",
        title = L"\text{Projekcja PCA } (%$subtitle)",
    )
    if n_show >= 2
        scatter!(
            ax_proj,
            transformed[:, 1],
            transformed[:, 2];
            markersize = 4.5,
            color = (palette[1], 0.75),
        )
    elseif n_show == 1
        scatter!(
            ax_proj,
            transformed[:, 1],
            zeros(size(transformed, 1));
            markersize = 4.5,
            color = (palette[1], 0.75),
        )
    end

    ax_evr = Axis(
        fig[1, 2],
        xlabel = L"\text{Główna składowa}",
        ylabel = L"\text{EVR}",
        limits = (0.5, max(length(evr), 1) + 0.5, 0, 1),
        title = L"\text{Współczynnik wariancji wyjaśnionej}",
    )
    if !isempty(evr)
        barplot!(ax_evr, 1:length(evr), evr; color = palette[2])
    end

    return fig
end

function plot_lle_dim(dataset::AbstractMatrix{<:Real}, k::Int, d::Int, tau::Real)
    set_publication_theme()

    palette = Makie.theme(:Palette, :color)[]
    lle_data = run_lle_per_time(dataset; k = k, d = d)

    if !haskey(lle_data.lle_results, tau)
        error("Wartość tau = $tau nie została znaleziona w wynikach LLE.")
    end

    embedding = lle_data.lle_results[tau]

    fig = Figure(size = (600, 500))
    ax = Axis(
        fig[1, 1], 
        title = L"\text{LLE: } k=%$k, d=%$d, \tau=%$tau"
    )

    if d == 1
        ax.xlabel = L"\text{LLE1}"
        ax.ylabel = L"\text{Wartość stała}"
        scatter!(
            ax,
            embedding[:, 1],
            zeros(size(embedding, 1));
            markersize = 4.5,
            color = (palette[1], 0.75),
        )
    else
        ax.xlabel = L"\text{LLE1}"
        ax.ylabel = L"\text{LLE2}"
        scatter!(
            ax,
            embedding[:, 1],
            embedding[:, 2];
            markersize = 4.5,
            color = (palette[1], 0.75),
        )
    end

    return fig
end

function plot_lle_dim!(ax::Axis, dataset::AbstractMatrix{<:Real}, k::Int, d::Int, tau::Real)
    set_publication_theme()

    palette = Makie.theme(:Palette, :color)[]
    lle_data = run_lle_per_time(dataset; k = k, d = d)
    embedding = lle_data.lle_results[tau]

    ax.title = L"\text{LLE: } k=%$k, d=%$d, \tau=%$tau"

    if d == 1
        ax.xlabel = L"\text{Odwzorowanie 1}"
        ax.ylabel = L"\text{Wartość stała}"
        scatter!(
            ax,
            embedding[:, 1],
            zeros(size(embedding, 1));
            markersize = 4.5,
            color = (palette[1], 0.75),
        )
    else
        ax.xlabel = L"\text{Odwzorowanie LLE1}"
        ax.ylabel = L"\text{Odwzorowanie LLE2}"
        scatter!(
            ax,
            embedding[:, 1],
            embedding[:, 2];
            markersize = 4.5,
            color = (palette[1], 0.75),
        )
    end
    return nothing
end

function plot_simulation_lle(dataset::AbstractMatrix{<:Real}, k::Int, d::Int, tau_zakres)
    set_publication_theme()

    palette = Makie.theme(:Palette, :color)[]
    liczba_wykresow = length(tau_zakres)
    kolumny = 2
    wiersze = ceil(Int, liczba_wykresow / kolumny)

    bok_kwadratu = max(kolumny, wiersze)

    fig = Figure(size = (bok_kwadratu * 400, bok_kwadratu * 400))

    for (i, tau) in enumerate(tau_zakres)
        row = (i - 1) ÷ kolumny + 1
        col = (i - 1) % kolumny + 1

        ax = Axis(fig[row, col], aspect = 1)
        plot_lle_dim!(ax, dataset, k, d, tau)
    end

    return fig
end

function animate_pca_evolution(
    dataset::AbstractMatrix{<:Real};
    filename::String = "pca_evolution.gif",
    fps::Int = 15,
    n_components::Int = 2,
    method::Symbol = :minmax,
    gamma::Float64 = 1.0,
    tau_tol::Float64 = 1e-8
)
    set_publication_theme()

    palette = Makie.theme(:Palette, :color)[]
    taus = sort(unique(dataset[:, 1]))

    fig = Figure(size = (800, 600))

    title_obs = Observable(L"\text{Projekcja PCA } (\tau = %$(taus[1]))")
    ax = Axis(
        fig[1, 1],
        xlabel = L"\text{PC1}",
        ylabel = L"\text{PC2}",
        title = title_obs
    )

    pts_obs = Observable(Point2f[])

    scatter!(
        ax,
        pts_obs;
        markersize = 6.0,
        color = (palette[1], 0.75)
    )

    record(fig, filename, taus; framerate = fps) do t
        val = round(t, digits=3)
        title_obs[] = L"\text{Projekcja PCA } (\tau = %$(val))"

        d = abs.(dataset[:, 1] .- t)
        mask = d .<= tau_tol
        data_slice = dataset[mask, :]

        features = data_slice[:, 2:3]

        if size(features, 1) > 1
            pca_result = if method === :minmax
                run_pca(features; n_components = n_components)
            elseif method === :kernel
                run_pca_kernel(features; n_components = n_components, gamma = gamma)
            else
                error("Unknown PCA method. Choose :minmax or :kernel.")
            end

            transformed = pca_result.transformed

            if size(transformed, 2) >= 2
                pts_obs[] = [Point2f(transformed[i, 1], transformed[i, 2]) for i in 1:size(transformed, 1)]
            else
                pts_obs[] = [Point2f(transformed[i, 1], 0.0) for i in 1:size(transformed, 1)]
            end

            autolimits!(ax)
        end
    end

    return filename
end
