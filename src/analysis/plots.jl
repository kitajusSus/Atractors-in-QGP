using ColorSchemes
using GLMakie
import GLMakie: Axis
using LaTeXStrings
# import Colors


"""
    set_publication_theme(; cmap = :devon, n_colors = 10, bg_color = RGBf(0.96, 0.96, 0.96))

`cmap` - np. :davos, :lajolla, :devon, :haline, :phase

"""
function set_publication_theme(;
        cmap = :devon,         # :davos, :lajolla, :devon, :haline, :phase
        n_colors = 25,         #
        bg_color = RGBf(0.98, 0.98, 0.98)
    )

    scientific_palette = Makie.resample_cmap(cmap, n_colors)

    return set_theme!(
        Theme(
            font = "Libertinus Serif",
            fontsize = 24,
            figure_padding = 20,

            Axis = (
                backgroundcolor = bg_color,
                titlesize = 28,
                xlabelsize = 26,
                ylabelsize = 26,
                xticklabelsize = 20,
                yticklabelsize = 20,

                xgridstyle = :dash,
                ygridstyle = :dash,
                xgridcolor = RGBAf(0.8, 0.8, 0.8, 0.7),
                ygridcolor = RGBAf(0.8, 0.8, 0.8, 0.7),

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

                # xminorticksvisible = true,
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
                backgroundcolor = RGBAf(1.0, 1.0, 1.0, 0.9),
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
                # strokewidth = 0.3,
                # strokecolor = :white,
            ),

            Palette = (
                color = scientific_palette,
                patchcolor = scientific_palette,
            ),
        )
    )
end


const PLOT_KEYS = Dict(
    :T => (L"T\,[\mathrm{fm}^{-1}]", (x, _) -> x[2]),
    :A => (L"\mathcal{A}", (x, _) -> x[3]),
    :tauT => (L"\tau T", (x, _) -> x[1] * x[2]),
    :Tdot => (L"\dot{T}\,[\mathrm{fm}^{-2}]", (x, _) -> (x[2] / x[1]) * (-1 / 3 + x[3] / 18)),
    :tau2Tdot => (L"\tau^2 \dot{T}", (x, _) -> x[1]^2 * ((x[2] / x[1]) * (-1 / 3 + x[3] / 18))),
    # testowanie normalizowanie
    :T_norm => (
        L"T / T_{\mathrm{max}}",
        (x, slice) -> x[2] / maximum(slice[:, 2]),
    ),
    :A_norm => (
        L"\mathcal{A} / \mathcal{A}_{\mathrm{max}}",
        (x, slice) -> x[3] / maximum(slice[:, 3]),
    ),
    # lambda funkcje kongo lekkie
    # Tdot / Tdot_max dla danej chwili czasu
    :Tdot_norm => (
        L"\dot{T} / \dot{T}_{\mathrm{max}}",
        (x, slice) -> begin
            Tdot_all = (slice[:, 2] ./ slice[:, 1]) .* (-1 / 3 .+ slice[:, 3] ./ 18)

            Tdot_current = (x[2] / x[1]) * (-1 / 3 + x[3] / 18)

            return Tdot_current / maximum(Tdot_all)
        end,
    ),
    # pod konkretne publikacje tutaj 2020 Hydrodynamics in Phase Space 0.22 to mój czas
    # inizjalitacji τ₀
    :tauT_2020 => (L"\tau_0 T", (x, _) -> 0.2 * x[2]),
    :tau2Tdot_2020 => (L"\tau_0^2 \dot{T}", (x, _) -> 0.2^2 * ((x[2] / x[1]) * (-1 / 3 + x[3] / 18)))
)

function resolve_def(def)
    if def isa Symbol
        @assert haskey(PLOT_KEYS, def) "Unknown plot key: $def"
        return PLOT_KEYS[def]
    end

    if isa(def, Tuple)&& length(def) == 2
        return (def[1], def[2])
    end

    error("Axis definition must be Symbol or Tuple(label, function).")
end

function get_data(dataset::AbstractMatrix{<:Real}, t::Real, xdef, ydef; is_attractor = false)
    xlbl, xfn = resolve_def(xdef)
    ylbl, yfn = resolve_def(ydef)

    if is_attractor
        selected = dataset
    else
        rows = findall(isapprox.(dataset[:, 1], t; atol = 1.0e-8))
        if isempty(rows)
            nearest = argmin(abs.(dataset[:, 1] .- t))
            rows = findall(isapprox.(dataset[:, 1], dataset[nearest, 1]; atol = 1.0e-8))
        end
        selected = dataset[rows, :]
    end

    x = [xfn(selected[i, :], selected) for i in 1:size(selected, 1)]
    y = [yfn(selected[i, :], selected) for i in 1:size(selected, 1)]

    return (x = x, y = y, xlabel = xlbl, ylabel = ylbl)
end


"""
    get_limits(dataset::AbstractMatrix{<:Real}, def; times = nothing, padding = 0.05)

Calculates the padded `(min, max)` range for a single plot key definition `def`.
If `times` is `nothing`, it uses all unique `tau` values in the dataset.
"""
function get_limits(dataset::AbstractMatrix{<:Real}, def; times = nothing, padding = 0.05)
    lbl, fn = resolve_def(def)
    ts = isnothing(times) ? unique(dataset[:, 1]) : times
    vals = Float64[]
    for t in ts
        rows = findall(isapprox.(dataset[:, 1], t; atol = 1.0e-8))
        if isempty(rows)
            nearest = argmin(abs.(dataset[:, 1] .- t))
            rows = findall(isapprox.(dataset[:, 1], dataset[nearest, 1]; atol = 1.0e-8))
        end
        selected = dataset[rows, :]
        for i in 1:size(selected, 1)
            push!(vals, fn(selected[i, :], selected))
        end
    end
    if isempty(vals)
        return (0.0, 1.0)
    end
    val_min, val_max = minimum(vals), maximum(vals)
    dval = val_max - val_min
    if dval ≈ 0.0
        dval = val_min ≈ 0.0 ? 1.0 : abs(val_min) * 0.1
    end
    return (val_min - padding * dval, val_max + padding * dval)
end

"""
    get_limits(dataset::AbstractMatrix{<:Real}, xdef, ydef; times = nothing, padding = 0.05)

Calculates the padded `(xmin, xmax, ymin, ymax)` ranges for both axis definitions `xdef` and `ydef`.
If `times` is `nothing`, it uses all unique `tau` values in the dataset.
"""
function get_limits(dataset::AbstractMatrix{<:Real}, xdef, ydef; times = nothing, padding = 0.05)
    xlbl, xfn = resolve_def(xdef)
    ylbl, yfn = resolve_def(ydef)
    ts = isnothing(times) ? unique(dataset[:, 1]) : times
    x_vals = Float64[]
    y_vals = Float64[]
    for t in ts
        rows = findall(isapprox.(dataset[:, 1], t; atol = 1.0e-8))
        if isempty(rows)
            nearest = argmin(abs.(dataset[:, 1] .- t))
            rows = findall(isapprox.(dataset[:, 1], dataset[nearest, 1]; atol = 1.0e-8))
        end
        selected = dataset[rows, :]
        for i in 1:size(selected, 1)
            row = selected[i, :]
            push!(x_vals, xfn(row, selected))
            push!(y_vals, yfn(row, selected))
        end
    end
    if isempty(x_vals) || isempty(y_vals)
        return (0.0, 1.0, 0.0, 1.0)
    end
    xmin, xmax = minimum(x_vals), maximum(x_vals)
    ymin, ymax = minimum(y_vals), maximum(y_vals)
    dx = xmax - xmin
    dy = ymax - ymin
    if dx ≈ 0.0
        dx = xmin ≈ 0.0 ? 1.0 : abs(xmin) * 0.1
    end
    if dy ≈ 0.0
        dy = ymin ≈ 0.0 ? 1.0 : abs(ymin) * 0.1
    end
    return (
        xmin - padding * dx, xmax + padding * dx,
        ymin - padding * dy, ymax + padding * dy,
    )
end

function _split_trajectories(dataset::AbstractMatrix{<:Real})
    @assert size(dataset, 2) >= 2 "Dataset must have at least columns [tau, feature1]."
    if size(dataset, 1) == 0
        return UnitRange{Int}[]
    end

    starts = Int[1]
    for i in 2:size(dataset, 1)
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

function plot_phase_space_grid(
        dataset::AbstractMatrix{<:Real},
        times,
        xdef,
        ydef;
        limits = nothing,
        attractor = nothing,
        attractor_points::Int = 150
    )
    set_publication_theme()
    palette = Makie.theme(:Palette).color[]

    n = length(times)
    ncols = min(3, n)
    nrows = ceil(Int, n / ncols)
    fig = Figure(size = (400 * ncols, 350 * nrows))

    lims = get_limits(dataset, xdef, ydef; times = times)

    for (i, t) in enumerate(times)
        row = (i - 1) ÷ ncols + 1
        col = (i - 1) % ncols + 1

        d = get_data(dataset, t, xdef, ydef)
        if i > 3

            ax = Axis(
                fig[row, col],
                title = L"\tau = %$(round(t, digits=2))\,\mathrm{fm}/c",
                xlabel = d.xlabel,
                ylabel = d.ylabel,
                limits = lims
            )
            ylims!(ax, 0, 4)
        else
            ax = Axis(
                fig[row, col],
                title = L"\tau = %$(round(t, digits=2))\,\mathrm{fm}/c",
                xlabel = d.xlabel,
                ylabel = d.ylabel,
                limits = lims
            )
        end

        # musi być najpiew by był pod kropkami
        # jeśli jest dany jeśnie nie to nie
        if !isnothing(attractor)
            attr_line = get_attractor_line_for_frame(
                dataset,
                attractor,
                t,
                xdef,
                ydef;
                limits = [0, lims[2], lims[3], lims[4]],
                n_points = attractor_points,
            )
            xlims!(ax, 0, lims[2])
            lines!(ax, attr_line.x, attr_line.y; color = (RGBAf(176 / 255, 0, 0, 1.0), 0.86), linewidth = 5, label = "Atraktor")
        end
        scatter!(
            ax, d.x[1:1000], d.y[1:1000]; markersize = 3.5, color = (palette[2], 0.9),
            strokecolor = (palette[6], 0.01),
            strokewidth = 0.1
        )
    end
    return fig
end

function plot_attractor(
        dataset::AbstractMatrix{<:Real},
        T0_target::Union{Real, AbstractVector{<:Real}, Nothing} = nothing;
        attractor::Union{Nothing, AbstractMatrix{<:Real}} = nothing,
        tol::Real = 0.05,
        n_trajectories::Int = 50,
        group_tol::Real = 0.01
    )
    set_publication_theme()
    fig = Figure(size = (950, 600))
    ax = Axis(
        fig[1, 1],
        xlabel = L"w = \tau T",
        ylabel = L" \mathcal{A}",
        limits = (0, 3, -1, 8)
    )

    trajs = _split_trajectories(dataset)

    if T0_target isa Real
        wybrane_trajs = [tr for tr in trajs if abs(dataset[tr[1], 2] - T0_target) < tol]
    elseif T0_target isa AbstractVector
        wybrane_trajs = [tr for tr in trajs if any(abs(dataset[tr[1], 2] - t) < tol for t in T0_target)]
    else
        wybrane_trajs = trajs
    end
    sort!(wybrane_trajs, by = tr -> dataset[tr[1], 2])

    if !isnothing(attractor)
        omega_attr = attractor[:, 1] .* attractor[:, 2]
        A_attr = attractor[:, 3]

        lines!(
            ax, omega_attr, A_attr,
            color = (:red, 1),
            linewidth = 7.0,
            label = L"\text{Teoretyczny atraktor}"
        )
        axislegend(ax, position = :rt)
    end
    step = max(1, length(wybrane_trajs) ÷ n_trajectories)
    selected_trajs = wybrane_trajs[1:step:end]

    T0s = [dataset[tr[1], 2] for tr in selected_trajs]
    min_T0, max_T0 = extrema(T0s)
    if min_T0 == max_T0
        min_T0 -= 0.1
        max_T0 += 0.1
    end

    for tr in selected_trajs
        T0 = dataset[tr[1], 2]
        norm_val = (T0 - min_T0) / (max_T0 - min_T0)
        col = get(ColorSchemes.batlow, norm_val)

        omega_traj = dataset[tr, 1] .* dataset[tr, 2]
        A_traj = dataset[tr, 3]

        lines!(
            ax, omega_traj, A_traj,
            color = (col, 0.5),
            linewidth = 1.5,
        )
    end


    hlines!(ax, [0.0], color = :black, linestyle = :dash, linewidth = 1.5)

    Colorbar(fig[1, 2], colormap = :batlow, limits = (min_T0, max_T0))

    return fig
end

function plot_thermodynamics_evolution(dataset::AbstractMatrix{<:Real})
    set_publication_theme()

    palette = Makie.theme(:Palette).color[]
    trajs = _split_trajectories(dataset)
    trajektorie = trajs[1:10:end]
    fig = Figure(size = (1200, 700))
    ax1 = Axis(
        fig[1, 1],
        title = L"\text{Ewolucja Temperatury } T\,[\mathrm{fm}^{-1}]\; \text{w czasie własnym } \tau",
        xlabel = L"\tau\,[\mathrm{fm}/c]",
        ylabel = L"T\,[\mathrm{fm}^{-1}]",
    )

    for tr in trajektorie
        lines!(
            ax1,
            dataset[tr, 1],
            dataset[tr, 2],
            color = (palette[3], 0.2),
            linewidth = 1.5,
        )
    end

    ax2 = Axis(
        fig[2, 1],
        title = L"\text{Ewolucja Anizotropii}\; \mathcal{A(τ)}\; \text{ w czasie własnym}",
        xlabel = L"\tau\,[\mathrm{fm}/c]",
        ylabel = L"\mathcal{A}",
    )

    for tr in trajektorie
        lines!(
            ax2,
            dataset[tr, 1],
            dataset[tr, 3],
            color = (palette[2], 0.2),
            linewidth = 1.5,
        )
    end

    hlines!(
        ax2,
        [0.0],
        color = :red,
        linestyle = :dash,
        linewidth = 2.0,
        label = L"\mathcal{A}=0\;(\text{Anizotropia} = 0)",
    )
    axislegend(ax2, position = :rt)

    linkxaxes!(ax1, ax2)
    return fig
end


function plot_phase_space_evolution(dataset::AbstractMatrix{<:Real})
    set_publication_theme()
    palette = Makie.theme(:Palette).color[]
    trajs = _split_trajectories(dataset)
    trajektorie = trajs[1:10:end]
    fig = Figure(size = (950, 620))
    ax = Axis(
        fig[1, 1],
        title = L"\text{Ewolucja w przestrzeni fazowej } (T, \mathcal{A})\; \text{w czasie własnym } \tau",
        xlabel = L"T\,[\mathrm{fm}^{-1}]",

        ylabel = L"\mathcal{A}",
    )
    for (i, tr) in enumerate(trajektorie)
        col = palette[mod1(i, length(palette))]
        lines!(
            ax,
            dataset[tr, 2],
            dataset[tr, 3],
            color = (col, 0.75),
            linewidth = 1.5,
        )
    end
    hlines!(
        ax,
        [0.0],
        color = :black,
        linestyle = :dash,
        linewidth = 2.0,
        label = L"\mathcal{A}=0\;(\text{Anizotropia} = 0)",
    )
    axislegend(ax, position = :rt)
    return fig
end


function plot_attractor_Aw_T(
        dataset::AbstractMatrix{<:Real},
        T_zadana::Real;
        attractor::Union{Nothing, AbstractMatrix{<:Real}} = nothing,
        T_tol::Real = 1.0e-3
    )
    set_publication_theme()

    T_col = dataset[:, 2]
    rows = findall(isapprox.(T_col, T_zadana; atol = T_tol))
    if isempty(rows)
        nearest_idx = argmin(abs.(T_col .- T_zadana))
        T_najblizsze = T_col[nearest_idx]
        rows = findall(isapprox.(T_col, T_najblizsze; atol = T_tol))
        println("Brak dokładnego T = $T_zadana w danych. Użyto najbliższego znalezionego: T ≈ $(round(T_najblizsze, digits = 4))")
    end
    slice = dataset[rows, :]

    omega_slice = slice[:, 1] .* slice[:, 2]
    A_slice = slice[:, 3]

    fig = Figure(size = (850, 600))
    ax = Axis(
        fig[1, 1],
        title = L"\text{Atraktor }\mathcal{A}(\omega)\text{ dla } T = %$(round(T_zadana, digits=3))\,\mathrm{fm}^{-1}",
        xlabel = L"\text{Czas uniwersalny } \omega = \tau T",
        ylabel = L"\text{Anizotropia } \mathcal{A} (w)"
    )

    if !isnothing(attractor)
        omega_attr = attractor[:, 1] .* attractor[:, 2]
        A_attr = attractor[:, 3]
        lines!(
            ax, omega_attr, A_attr,
            color = :red,
            linewidth = 3.0,
            label = L"\text{Teoretyczny atraktor}"
        )
    end

    # Rysujemy chmurę punktów z symulacji
    palette = Makie.theme(:Palette).color[]
    scatter!(
        ax, omega_slice, A_slice,
        markersize = 8,
        color = (palette[1], 0.6),
        strokecolor = :white,
        strokewidth = 0.5,
        label = L"\text{Dane symulacji}"
    )

    axislegend(ax, position = :rb)

    return fig
end


function plot_phase_space_evolution_3d(dataset::AbstractMatrix{<:Real})
    set_publication_theme()
    palette = Makie.theme(:Palette).color[]
    trajs = _split_trajectories(dataset)
    trajektorie = trajs[1:10:end]
    fig = Figure(size = (950, 750))
    ax = Axis3(
        fig[1, 1],
        title = L"\text{Ewolucja w przestrzeni fazowej } (T, \mathcal{A}, \tau)",
        xlabel = L"T\,[\mathrm{fm}^{-1}]",
        ylabel = L"\mathcal{A}",
        zlabel = L"\tau\,[\mathrm{fm}]",
        azimuth = 1.3π,
        elevation = 0.15π,
    )
    for (i, tr) in enumerate(trajektorie)
        col = palette[mod1(i, length(palette))]
        lines!(
            ax,
            dataset[tr, 2],
            dataset[tr, 3],
            dataset[tr, 1],
            color = (col, 0.7),
            linewidth = 1.5,
        )
    end
    T_range = range(
        minimum(dataset[:, 2]),
        maximum(dataset[:, 2]),
        length = 2,
    )
    τ_range = range(
        minimum(dataset[:, 1]),
        maximum(dataset[:, 1]),
        length = 2,
    )
    surface!(
        ax,
        T_range,
        zeros(2, 2),         # A=0
        repeat(τ_range, 1, 2)',
        color = fill((:gray, 0.15), 2, 2),
        transparency = true,
    )
    return fig
end


function plot_pca_evr_over_time(
        dataset::AbstractMatrix{<:Real};
        n_components::Int = 2,
        method::Symbol = :minmax,
        gamma::Float64 = 1.0,
        feature_cols::AbstractVector = collect(2:size(dataset, 2)),
        plot_title::Union{String, LaTeXString} = L"\text{Explained Variance Ratio (EVR) w funkcji czasu}",
        x_label::Union{String, LaTeXString} = L"\tau\,[\mathrm{fm}/c]",
        y_label::Union{String, LaTeXString} = L"\text{EVR}",
        tau_min::Union{Real, Nothing} = nothing
    )
    set_publication_theme()

    result = run_pca_per_time(
        dataset;
        n_components = n_components,
        method = method,
        gamma = gamma,
        feature_cols = feature_cols,
    )
    taus = result.taus
    evr = result.explained_variance_ratio

    t_min = isnothing(tau_min) ? minimum(taus) : Float64(tau_min)

    fig = Figure(size = (850, 600))
    ax = Axis(
        fig[1, 1],
        title = plot_title,
        xlabel = x_label,
        ylabel = y_label,

        limits = (t_min, maximum(taus), 0, 1.05),
    )

    hlines!(ax, [1.0], color = :gray45, linestyle = :dash, label = L"100\%")

    palette = Makie.theme(:Palette).color[]
    kolory = [palette[3], palette[7], palette[3], palette[7], palette[2], palette[8]]
    for comp in 1:n_components
        vals = evr[:, comp]
        mask = .!isnan.(vals)
        if any(mask)
            lines!(
                ax,
                taus[mask],
                vals[mask],
                linewidth = 3.0,
                color = kolory[mod1(comp, length(kolory))],
                label = L"\text{PC}%$(comp)",
            )
            # if comp == 1
            #     band!(
            #         ax,
            #         taus[mask],
            #         zeros(sum(mask)),
            #         vals[mask],
            #         color = (palette[1], 0.15),
            #     )
            # end
        end
    end

    axislegend(ax, position = :rb)
    return fig
end

function plot_pca_bar_variance(
        dataset::AbstractMatrix{<:Real};
        tau::Real = 1.14,
        method::Symbol = :minmax,
        gamma::Float64 = 1.0
    )
    set_publication_theme()
    palette = Makie.theme(:Palette).color[]

    Xtau = if size(dataset, 2) == 3
        _, sliced = get_tau_slice(dataset, tau; atol = 1.0e-8, feature_cols = [2, 3])
        sliced
    else
        Matrix{Float64}(dataset)
    end

    pca_result = if method === :minmax
        run_pca(Xtau; n_components = 2)
    elseif method === :kernel
        run_pca_kernel(Xtau; n_components = 2, gamma = gamma)
    else
        error("Nieznana metoda PCA. Wybierz :minmax lub :kernel.")
    end

    evr = pca_result.explained_variance_ratio_full
    n_comp = length(evr)
    cumulative_evr = cumsum(evr)

    fig = Figure(size = (850, 550))
    ax = Axis(
        fig[1, 1],
        title = L"\text{Analiza Wariancji — — PCA z normalizacją min-max}",
        xlabel = L"\text{Główna Składowa}",
        ylabel = L"\text{EVR [\%]}",
        xticks = (1:n_comp, [L"\mathrm{PC}_%$(i)" for i in 1:n_comp]),
        limits = (0.4, n_comp + 0.6, 0, 110)
    )

    evr_percent = evr .* 100
    barplot!(ax, 1:n_comp, evr_percent, color = palette[1], width = 0.5)

    for i in 1:n_comp
        text!(
            ax,
            i,
            evr_percent[i] + 3,
            text = string(round(evr_percent[i], digits = 1), "%"),
            align = (:center, :bottom),
            fontsize = 18,
            font = "Libertinus Serif"
        )
    end

    cum_percent = cumulative_evr .* 100
    lines!(ax, 1:n_comp, cum_percent, color = :gray35, linestyle = :dash, linewidth = 2.5)
    scatter!(ax, 1:n_comp, cum_percent, color = palette[2], markersize = 12, label = L"\text{Suma EVR}")

    axislegend(ax, position = :rc)

    return fig
end

function plot_pca_summary(
        dataset::AbstractMatrix{<:Real};
        tau::Union{Nothing, Real} = nothing,
        tau_tol::Float64 = 1.0e-8,
        tau_mode::Symbol = :nearest,
        n_components::Int = 2,
        method::Symbol = :minmax,
        gamma::Float64 = 1.0,
    )

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
function animate_pca_evolution(
        dataset::AbstractMatrix{<:Real};
        filename::String = "pca_evolution.gif",
        fps::Int = 15,
        n_components::Int = 2,
        method::Symbol = :minmax,
        gamma::Float64 = 1.0,
        tau_tol::Float64 = 1.0e-8
    )
    set_publication_theme()

    palette = Makie.theme(:Palette).color[]
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
        val = round(t, digits = 3)
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

            reset_limits!(ax)
        end
    end

    return filename
end
###########################
### FUNKCJE DO LLE
###############################

"""
Analiza spektrum wartości własnych dla danego czasu 
    
    function plot_lle_spectrum_statistics(dataset; τ, k_values = 5:5:50, ile_λ = 4)

- `ile_λ` - ile pierwszych wartości własnych pokazać na wykresie (domyślnie 4) 

ale trzeba pamiętać że pierwsza wartość własna jest (i powinna być) zawsze λ₁ = 0

"""
function plot_lle_spectrum_statistics(dataset; τ, k_values = 5:5:50, ile_λ = 4)
    set_publication_theme()

    palette = Makie.theme(:Palette).color[]
    spectra, ks = lle_spectrum_over_k(dataset; tau = τ, k_values = k_values)
    μ, σ = spectrum_statistics(spectra)

    x = 1:ile_λ
    k_min = minimum(k_values)
    k_max = maximum(k_values)

    fig = Figure(size = (1400, 700))
    ax = Axis(
        fig[1, 1],
        title = L"\text{LLE analiza wartości własnych [K z zakresu (%$k_min - %$k_max)] (średnia } \pm \sigma \text{) } \tau=%$τ",
        xlabel = L"\text{Indeks wartości własnej } \lambda_{i}",
        ylabel = L"\text{Wartość } \lambda_{i}"
    )

    scatterlines!(
        ax, x, μ[1:ile_λ],
        linewidth = 2,
        color = (palette[5], 0.3),
        markercolor = palette[5],
        markersize = 15,
        label = L"\text{Średnia } i\text{-ta } \lambda_{i}"
    )

    band!(
        ax,
        x,
        μ[1:ile_λ] .- σ[1:ile_λ],
        μ[1:ile_λ] .+ σ[1:ile_λ],
        color = (palette[3], 0.2)
    )

    axislegend(ax, position = :lt)

    return fig
end


"""
    plot_lle_spectrum_statistics_grid(dataset, taus; k_values = 5:5:50, ile_λ = 4, atol = 1.0e-8)
"""
function plot_lle_spectrum_statistics_grid(dataset, taus; k_values = 5:5:50, ile_λ = 4, atol = 1.0e-8)
    set_publication_theme()
    palette = Makie.theme(:Palette).color[]

    n = length(taus)
    ncols = min(3, n)
    nrows = ceil(Int, n / ncols)
    spectra1, ks1 = lle_spectrum_over_k(dataset; tau = taus[1], k_values = k_values, atol = atol)
    μ1, σ1 = spectrum_statistics(spectra1)

    y_min = minimum(μ1[1:ile_λ] .- σ1[1:ile_λ])
    y_max = maximum(μ1[1:ile_λ] .+ σ1[1:ile_λ])
    y_pad = (y_max - y_min) * 0.07
    lims = (0.9, ile_λ + 0.05, y_min - y_pad, y_max + y_pad)

    fig = Figure(size = (450 * ncols + 40, 350 * nrows + 120))

    k_min = minimum(k_values)
    k_max = maximum(k_values)
    title_text = L"\text{LLE analiza wartości własnych [K z zakresu (%$k_min - %$k_max)] (średnia } \pm \sigma \text{)}"
    Label(fig[1, 2:(ncols + 1)], title_text, fontsize = 36, font = :bold, padding = (0, 0, 10, 0))

    Label(fig[2:(nrows + 1), 1], L"\text{Wartość } \lambda_{i}", rotation = pi / 2, font = :bold, fontsize = 32)
    Label(fig[nrows + 2, 2:(ncols + 1)], L"\text{Indeks wartości własnej } \lambda_{i}", font = :bold, fontsize = 32)

    x = 1:ile_λ

    for (i, tau) in enumerate(taus)
        row = (i - 1) ÷ ncols + 2
        col = (i - 1) % ncols + 2
        if i == 1
            μ, σ = μ1, σ1
        else
            spectra, ks = lle_spectrum_over_k(dataset; tau = tau, k_values = k_values, atol = atol)
            μ, σ = spectrum_statistics(spectra)
        end

        ax = Axis(
            fig[row, col],
            title = L"\tau = %$tau\,\mathrm{fm}/c",
            limits = lims,
            xticks = 1:ile_λ
        )

        scatterlines!(
            ax, x, μ[1:ile_λ],
            linewidth = 2,
            color = (palette[5], 0.3),
            markercolor = palette[5],
            markersize = 15,
            label = L"\text{Średnia } i\text{-ta } \lambda_{i}"
        )

        band!(
            ax,
            x,
            μ[1:ile_λ] .- σ[1:ile_λ],
            μ[1:ile_λ] .+ σ[1:ile_λ],
            color = (palette[3], 0.2)
        )

        axislegend(ax, position = :lt)

        row_idx = (i - 1) ÷ ncols + 1
        col_idx = (i - 1) % ncols + 1
        is_bottom = (row_idx == nrows) || (i + ncols > n)
        is_left = (col_idx == 1)

        if !is_bottom
            hidexdecorations!(ax, grid = false)
        end
        if !is_left
            hideydecorations!(ax, grid = false)
        end
        println("Wyliczone wartości własne $tau: ", μ[1:ile_λ])
    end

    return fig
end

@views function plot_lle_σ_spectrum(dataset, taus; k_values = 5:5:50, ile_λ = 3, atol = 1.0e-8)
    set_publication_theme()
    palette = Makie.theme(:Palette).color[]
    n_taus = length(taus)
    sigma_matrix = fill(NaN, n_taus, ile_λ)

    for (i, tau) in enumerate(taus)
        spectra, ks = lle_spectrum_over_k(dataset; tau = tau, k_values = k_values, atol = atol)

        μ, σ = spectrum_statistics(spectra)
        N = min(ile_λ, length(σ))
        sigma_matrix[i, 1:N] .= σ[1:N] .* 10^(8)

    end

    fig = Figure(size = (800, 600))
    ax = Axis(
        fig[1, 1],
        title = L"\text{Zmiana } \sigma_i \text{ w funkcji } \tau \text{ dla } i \le %$ile_λ",
        xlabel = L"\tau \ [\mathrm{fm/c}]",
        ylabel = L"\sigma_i \cdot 10^{8}",

        xticks = 0.2:0.05:maximum(taus)
    )
    for j in 1:ile_λ
        xs = taus
        ys = sigma_matrix[:, j]
        lines!(ax, xs, ys, linewidth = 3, color = palette[(j + 1) * 3], label = L"\sigma_{%$j}")
        scatter!(ax, xs, ys, color = palette[j + 3], markersize = 5)
    end

    axislegend(ax, position = :rt)

    return fig
end

function plot_lle_dim(dataset::AbstractMatrix{<:Real}, k::Int, d::Int, tau::Real)
    set_publication_theme()

    palette = Makie.theme(:Palette, :color)[]
    lle_data = run_lle_per_time(dataset; k = k, d = d)
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

    palette = Makie.theme(:Palette).color[]
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

    palette = Makie.theme(:Palette).color[]
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

function plot_lle_embedding(dataset::AbstractMatrix{<:Real}; labels = nothing, title = L"\text{Projekcja LLE}", osie = nothing)
    set_publication_theme()
    palette = Makie.theme(:Palette).color[]
    dim = size(dataset, 2)

    fig = Figure(size = (750, 600))

    color_param = isnothing(labels) ? (palette[1], 0.75) : labels
    cmap_param = isnothing(labels) ? :viridis : :jet

    if dim >= 3
        ax = Axis3(fig[1, 1], title = title)
        #xlabel = L"\text{LLE1}", ylabel = L"\text{LLE2}", zlabel = L"\text{LLE3}")
        scatter!(ax, dataset[:, 1], dataset[:, 2], dataset[:, 3], color = color_param, colormap = cmap_param, markersize = 6)
    elseif dim == 2
        ax = Axis(fig[1, 1], title = title, xlabel = L"\text{LLE1}", ylabel = L"\text{LLE2}")
        scatter!(ax, dataset[:, 1], dataset[:, 2], color = color_param, colormap = cmap_param, markersize = 6)
    else
        ax = Axis(fig[1, 1], title = title, xlabel = L"\text{LLE1}")
        scatter!(ax, dataset[:, 1], zeros(size(dataset, 1)), color = color_param, colormap = cmap_param, markersize = 6)
    end

    return fig
end

function plot_lle_grid(lle_results::Dict, taus; labels = nothing)
    dim = size(first(values(lle_results)), 2)
    if dim >= 3
        return plot_lle_grid_3d(lle_results, taus; labels = labels)
    else
        return plot_lle_grid_2d(lle_results, taus; labels = labels)
    end
end

function plot_lle_grid_2d(lle_results::Dict, taus; labels = nothing)
    set_publication_theme()
    palette = Makie.theme(:Palette).color[]
    color_param = isnothing(labels) ? (palette[1], 0.75) : labels
    cmap_param = isnothing(labels) ? :viridis : :jet

    n = length(taus)
    ncols = min(3, n)
    nrows = ceil(Int, n / ncols)
    fig = Figure(size = (400 * ncols, 380 * nrows))

    for (i, t) in enumerate(taus)
        row = (i - 1) ÷ ncols + 1
        col = (i - 1) % ncols + 1
        embedding = lle_results[t]

        ax = Axis(fig[row, col], title = L"\tau = %$(round(t, digits=2))", xlabel = L"\text{LLE1}", ylabel = L"\text{LLE2}")
        scatter!(ax, embedding[:, 1], embedding[:, 2], color = color_param, colormap = cmap_param, markersize = 5)
    end

    return fig
end

function plot_lle_grid_3d(lle_results::Dict, taus; labels = nothing)
    set_publication_theme()
    palette = Makie.theme(:Palette).color[]
    color_param = isnothing(labels) ? (palette[1], 0.75) : labels
    cmap_param = isnothing(labels) ? :viridis : :jet

    n = length(taus)
    ncols = min(3, n)
    nrows = ceil(Int, n / ncols)
    fig = Figure(size = (400 * ncols, 380 * nrows))

    for (i, t) in enumerate(taus)
        row = (i - 1) ÷ ncols + 1
        col = (i - 1) % ncols + 1
        embedding = lle_results[t]

        ax = Axis3(fig[row, col], title = L"\tau = %$(round(t, digits=2))", xlabel = L"\text{LLE1}", ylabel = L"\text{LLE2}", zlabel = L"\text{LLE3}")
        scatter!(ax, embedding[:, 1], embedding[:, 2], embedding[:, 3], color = color_param, colormap = cmap_param, markersize = 5)
    end

    return fig
end


## funkcje do PINN
function plot_pinn_deff_evolution(results::AbstractMatrix{<:Real})
    set_publication_theme()
    palette = Makie.theme(:Palette).color[]
    fig = Figure(size = (850, 600))
    ax = Axis(
        fig[1, 1],
        title = L"\text{Ewolucja wymiaru efektywnego } d_{\mathrm{eff}} \text{ w czasie } \tau",
        xlabel = L"\tau\,[\mathrm{fm}/c]",
        ylabel = L"d_{\mathrm{eff}}",
        limits = (minimum(results[:, 1]), maximum(results[:, 1]), 0.95, maximum(results[:, 2]) * 1.05)
    )
    lines!(ax, results[:, 1], results[:, 2], color = palette[1], linewidth = 3.0)
    return fig
end


#############
## FUNKCJE TO NN
# ###
#
#
#
#
#


#########################
# funkcje do lid mle twonn, itd
# #########################
#
#
#


function plot_lid_dimension(
        dataset::AbstractMatrix{<:Real},
        #k::Int
    )
    set_publication_theme()

    palette = Makie.theme(:Palette).color[]
    fig = Figure(size = (1200, 700))
    ax = Axis(
        fig[1, 1],
        title = L"\text{Ewolucja estymowanego wymiaru LID w funkcji czasu dla różnych } k",
        xlabel = L"\text{Czas } \tau\,[\mathrm{fm}/c]",
        ylabel = L"\text{estymowany wymiar } d_{\mathrm{LID}}",
    )
    k_range = 1:10:100
    for (idx, current_k) in enumerate(k_range)
        taus, lid, _, _ = scan_intrinsic_dimensions(dataset, k = current_k)

        lines!(
            ax, taus, lid,
            color = palette[idx],
            linewidth = 2.5,
            label = latexstring("k = ", current_k)
        )

        # println("τ = $(taus[idx]), LID = $(lid[idx])")
        # println("Obliczony LID dla k = $current_k: ", lid)
    end
    return fig
end


function plot_twonn(dataset::AbstractMatrix{<:Real})
    set_publication_theme()

    palette = Makie.theme(:Palette).color[]
    fig = Figure(size = (1200, 700))
    ax = Axis(
        fig[1, 1],
        title = L"\text{Przykładowy wykres dla TWONN - Two Nearest Neighbors}",
        xlabel = L"\text{Czas } \tau\,[\mathrm{fm}/c]",
        ylabel = L"\text{estymowany wymiar } d_{\mathrm{TWONN}}",
    )
    taus, _, twonn, _ = scan_intrinsic_dimensions(dataset)
    lines!(ax, taus, twonn, color = palette[3], linewidth = 2.5, label = L"\text{TWONN}")
    axislegend(ax, position = :rt)
    return fig
end

function plot_lle_results_for_taus(
        dataset::AbstractMatrix{<:Real},
        target_taus::Vector{Float64};
        k::Integer = 20,
        d::Integer = 2,
        atol::Real = 1.0e-3,
        feature_cols::AbstractVector{<:Integer} = collect(2:size(dataset, 2))
    )

    set_publication_theme()
    palette = Makie.theme(:Palette).color[]

    results = run_lle_for_selected_taus(dataset, target_taus; k = k, d = d, atol = atol, feature_cols = feature_cols)
    lle_res = results.lle_results

    n = length(target_taus)
    ncols = min(3, n)
    nrows = ceil(Int, n / ncols)

    fig_width = max(1400, 450 * ncols)
    fig_height = max(700, 450 * nrows)
    fig = Figure(size = (fig_width, fig_height))

    for (i, t) in enumerate(target_taus)
        if !haskey(lle_res, t)
            continue
        end
        row = (i - 1) ÷ ncols + 1
        col = (i - 1) % ncols + 1

        embedding = lle_res[t]
        dim = size(embedding, 2)

        ax_title = L"\tau = %$(round(t, digits=2))"
        if dim >= 3
            ax = Axis3(fig[row, col], title = ax_title, xlabel = L"\text{LLE1}", ylabel = L"\text{LLE2}", zlabel = L"\text{LLE3}")
            scatter!(ax, embedding[:, 1], embedding[:, 2], embedding[:, 3], color = (palette[1], 0.75), markersize = 6)
        elseif dim == 2
            ax = Axis(fig[row, col], title = ax_title, xlabel = L"\text{LLE1}", ylabel = L"\text{LLE2}")
            scatter!(ax, embedding[:, 1], embedding[:, 2], color = (palette[1], 0.75), markersize = 6)
        else
            ax = Axis(fig[row, col], title = ax_title, xlabel = L"\text{LLE1}")
            scatter!(ax, embedding[:, 1], zeros(size(embedding, 1)), color = (palette[1], 0.75), markersize = 6)
        end
    end

    return fig
end

function plot_lle_spectrum_analysis(
        dataset::AbstractMatrix{<:Real},
        tau::Real;
        k_values = 5:5:50,
        atol::Real = 1.0e-8
    )

    set_publication_theme()
    palette = Makie.theme(:Palette).color[]

    spectra, ks = lle_spectrum_over_k(dataset; tau = tau, k_values = k_values, atol = atol)

    fig = Figure(size = (1400, 700))
    ax = Axis(
        fig[1, 1],
        title = L"\text{Widmo LLE dla różnych } k \text{, } \tau = %$(tau)",
        xlabel = L"\text{Indeks wartości własnej}",
        ylabel = L"\lambda"
    )

    for (i, (s, k)) in enumerate(zip(spectra, ks))
        col = palette[mod1(i, length(palette))]
        lines!(ax, 1:length(s), s, label = L"k = %$k", linewidth = 2.5, color = col)
    end

    axislegend(ax, position = :lt, nbanks = 2)
    return fig
end

function plot_lle_spectrum_scan_analysis(
        dataset::AbstractMatrix{<:Real},
        taus::Vector{Float64};
        k_values = 5:5:50
    )

    set_publication_theme()
    palette = Makie.theme(:Palette).color[]

    results = scan_lle_spectrum(dataset; taus = taus, k_values = k_values)

    fig = Figure(size = (1400, 700))
    ax = Axis(
        fig[1, 1],
        title = L"\text{Statystyki widma LLE w czasie}",
        xlabel = L"\text{Indeks wartości własnej}",
        ylabel = L"\lambda"
    )

    for (i, t) in enumerate(taus)
        if !haskey(results, t)
            continue
        end
        res = results[t]
        μ = res.mean
        σ = res.std
        x = 1:length(μ)

        col = palette[mod1(i, length(palette))]
        lines!(ax, x, μ, label = L"\tau = %$(t)", linewidth = 3, color = col)
        band!(ax, x, μ .- σ, μ .+ σ, color = (col, 0.2))
    end

    axislegend(ax, position = :rt, nbanks = 2)
    return fig
end


@views function plot_local_pca(
        dataset_loaded::AbstractMatrix{<:Real};
        n_slices::Int = 15,
        tablica_k::Vector{Int} = [10, 20, 40, 80, 160],
        feature_cols::AbstractVector{<:Integer} = 2:size(dataset_loaded, 2)
    )
    set_publication_theme()

    fig = Figure(size = (950, 600))
    ax = Axis(
        fig[1, 1],
        # title = L"\text{Wpływ parametru } K \text{ na działanie Algorytmu L-PCA jako } f(\tau)",
        xlabel = L"\tau\,[\mathrm{fm}/c]",
        ylabel = L"\text{Średni wymiar lokalny}",
        xticks = 0:0.25:maximum(dataset_loaded[:, 1]),
        yticks = 1:0.1:2,

        xautolimitmargin = (0.0, 0.05),
        yautolimitmargin = (0.05, 0.05)
    )

    palette = [:crimson, :dodgerblue, :forestgreen, :darkorange, :purple, :goldenrod, :darkcyan, :mediumvioletred]
    for (i, k_bazowe) in enumerate(tablica_k)

        k_drugie = k_bazowe .* 2
        tau_vals, mean_k1, std_k1 = compute_lpca(dataset_loaded, k_bazowe, n_slices, feature_cols = feature_cols)

        tau_drugie, mean_k2, std_k2 = compute_lpca(dataset_loaded, k_drugie, n_slices, feature_cols = feature_cols)
        c = palette[mod1(i, length(palette))]
        band!(
            ax,
            tau_vals,
            mean_k1,
            mean_k2;
            color = (c, 0.25)
        )

        lines!(
            ax,
            tau_vals,
            mean_k1;
            linewidth = 2.5,
            color = c,
            label = L"K = %$(k_bazowe)"
        )
    end

    axislegend(ax, position = :rt)

    return fig
end

function plot_map_lpca(
        dataset::AbstractMatrix{<:Real};
        zakres_K::AbstractVector{<:Integer} = [10, 20, 40, 80],
        n_slices::Int = 15,
        feature_cols::AbstractVector{<:Integer} = collect(2:size(dataset, 2))
    )

    set_publication_theme()
    palette = Makie.theme(:Palette).color[]

    results = [compute_lpca(dataset, k, n_slices; feature_cols) for k in zakres_K]
    tau_vals = results[1][1]
    dim_matrix = reduce(hcat, [r[2] for r in results])

    fig = Figure(size = (1400, 700))
    y_indices = [5:20:length(zakres_K)...]
    y_labels = string.(zakres_K[y_indices])

    ax = Axis(
        fig[1, 1],
        # title = L"\text{Mapa lokalnych wymiarów w funkcji } \mathcal{L_D}(\tau, K)",
        xlabel = L"\tau\,[\mathrm{fm}/c]",
        ylabel = L"K",
        xticks = 0:0.2:maximum(tau_vals),
        yticks = (y_indices, y_labels),
    )

    hm = heatmap!(
        ax, tau_vals, eachindex(zakres_K), dim_matrix;
        colormap = palette[8:end],
    )

    Colorbar(fig[1, 2], hm; label = L"\text{Lokalny wymiar}")

    return fig
end
