using Interpolations

using CSV
using DataFrames
using HDF5
using Serialization


function attractor_data(model::AbstractHydroModel; tau_max = 5.0, krok = 0.01, temperature_unit::Symbol = :fm)
    tspan_attr = (0.1, tau_max)

    ic_attr = [4, 20]
    sol_attr = solve_hydro(model, ic_attr, tspan_attr; saveat = krok)
    attractor_matrix = build_dataset([sol_attr]; temperature_unit = temperature_unit)

    return attractor_matrix
end


function attractor_data(
        model::MISModel;
        tau_max = 5.0,
        krok = 0.01,
        n_interp = 5000,
        w_min::Real = 1.0e-15
    )

    w_max = Float64(tau_max)
    p = model.params

    A0 = 6 * sqrt(p.eta_over_s / p.tau_pi)

    rhs_A(u, _, w) = begin
        A = u[1]
        dA = 18 * (8 * p.eta_over_s - w * A - (2 / 9) * p.tau_pi * A^2) /
            (p.tau_pi * w * (A + 12))
        return SVector(dA)
    end

    problem = ODEProblem(
        rhs_A,
        SVector(Float64(A0)),
        (Float64(w_min), w_max),
        nothing
    )

    sol = solve(problem, Rodas5(); abstol = 1.0e-9, reltol = 1.0e-9)

    ws = range(w_min, w_max, length = n_interp)

    A = [sol(w)[1] for w in ws]

    return hcat(collect(ws), ones(length(ws)), A)
end


"""
    get_attractor_line(...)

Project attractor onto arbitrary observables.
"""
function get_attractor_line(
        dataset::AbstractMatrix{<:Real},
        attractor::AbstractMatrix{<:Real},
        τ::Real,
        xdef,
        ydef;
        n_points::Int = 200
    )

    _, xmap = resolve_def(xdef)
    _, ymap = resolve_def(ydef)

    _, slice = get_tau_slice(
        dataset,
        τ;
        feature_cols = 1:size(dataset, 2)
    )

    Tmin = minimum(slice[:, 2])
    Tmax = maximum(slice[:, 2])

    Tgrid = range(
        Tmin * 0.0001,
        Tmax * 10;
        length = n_points
    )

    Aω = build_attractor_interpolant(attractor)

    x = Vector{Float64}(undef, n_points)
    y = Vector{Float64}(undef, n_points)

    ncols = size(dataset, 2)

    for (i, T) in enumerate(Tgrid)

        # T_fm = T > 50 ? T * FM_PER_MEV : T
        T_fm = T

        ω = τ * T_fm

        state = zeros(Float64, ncols)

        state[1] = τ
        state[2] = T
        state[3] = Aω(ω)

        x[i] = xmap(state, slice)
        y[i] = ymap(state, slice)

    end

    return (; x, y)
end

"""
    build_attractor_interpolant(attractor)

Construct an interpolation A(ω) from attractor data.
"""
function build_attractor_interpolant(
        attractor::AbstractMatrix{<:Real}
    )

    # T_is_mev = mean(attractor[:, 2]) > 50.0

    # Tfm = T_is_mev ?
    #     attractor[:, 2] .* FM_PER_MEV :
    #     attractor[:, 2]
    Tfm = attractor[:, 2]

    ω = attractor[:, 1] .* Tfm
    A = attractor[:, 3]

    p = sortperm(ω)

    ωsorted = ω[p]
    Asorted = A[p]

    return linear_interpolation(
        ωsorted,
        Asorted;
        extrapolation_bc = Line()
    )
end

"""
    temperature_grid(T_values::AbstractVector{<:Real}, n::Integer)

Construct a uniform temperature grid of length `n` from 95% of the minimum
to 105% of the maximum of `T_values`.
"""
function temperature_grid(T_values::AbstractVector{<:Real}, n::Integer)
    Tmin = minimum(T_values) * 0.95
    Tmax = maximum(T_values) * 1.05
    return temperature_grid(Tmin, Tmax, n)
end

function temperature_grid(Tmin::Real, Tmax::Real, n::Integer)
    @assert n > 0 "temperature grid must contain at least one point."
    if n == 1
        return [(Tmin + Tmax) / 2]
    end

    return range(
        Tmin,
        Tmax;
        length = n
    )
end

function interpolate_attractor_state(
        dataset::AbstractMatrix{S},
        attractor::AbstractMatrix{<:Real},
        τ::Real,
        T::Real
    ) where {S}
    omega = attractor[:, 1] .* attractor[:, 2]
    order = sortperm(omega)
    omega_sorted = omega[order]
    attr_sorted = attractor[order, :]

    target = τ * T
    state = Vector{S}(undef, size(dataset, 2))
    state[1] = τ
    state[2] = T

    n_copy = min(size(dataset, 2), size(attractor, 2))
    if target <= omega_sorted[1]
        state[3:n_copy] .= attr_sorted[1, 3:n_copy]
    elseif target >= omega_sorted[end]
        state[3:n_copy] .= attr_sorted[end, 3:n_copy]
    else
        hi = searchsortedfirst(omega_sorted, target)
        lo = hi - 1
        weight = (target - omega_sorted[lo]) / (omega_sorted[hi] - omega_sorted[lo])
        state[3:n_copy] .= (1 - weight) .* attr_sorted[lo, 3:n_copy] .+
            weight .* attr_sorted[hi, 3:n_copy]
    end

    if size(dataset, 2) > n_copy
        state[(n_copy + 1):end] .= zero(S)
    end
    return state
end

"""
    get_attractor_line_for_frame(dataset::AbstractMatrix{<:Real}, attractor_universe::AbstractMatrix{<:Real}, t::Real, xdef, ydef; limits=nothing, n_points=150)

Oblicza dokładne współrzędne linii teoretycznego atraktora dla wybranej klatki czasowej `t`,
dopasowując zakres i ewentualną normalizację do punktów numerycznych z symulacji (`dataset`).
"""
function get_attractor_line_for_frame(
        dataset::AbstractMatrix{<:Real},
        attractor_universe::AbstractMatrix{<:Real},
        t::Real,
        xdef,
        ydef;
        limits = nothing,
        n_points::Int = 150
    )
    _, xfn = resolve_def(xdef)
    _, yfn = resolve_def(ydef)

    _, cloud_slice = get_tau_slice(
        dataset,
        t;
        atol = 1.0e-8,
        feature_cols = 1:size(dataset, 2),
        nearest = true,
    )
    τ = cloud_slice[1, 1]

    # Helper function to convert axis value back to temperature T
    get_T_from_val = (val, def) -> begin
        if def === :T
            return val
        elseif def === :tauT
            return val / τ
        elseif def === :T_norm
            return val * maximum(cloud_slice[:, 2])
        elseif def === :tauT_heller
            return val / 0.2
        else
            return val
        end
    end

    Tmin, Tmax = if !isnothing(limits) && length(limits) == 4
        xmin, xmax, ymin, ymax = limits
        if xdef in (:T, :tauT, :T_norm, :tauT_heller)
            get_T_from_val(xmin, xdef), get_T_from_val(xmax, xdef)
        elseif ydef in (:T, :tauT, :T_norm, :tauT_heller)
            get_T_from_val(ymin, ydef), get_T_from_val(ymax, ydef)
        else
            minimum(dataset[:, 2]) * 0.95, maximum(dataset[:, 2]) * 1.05
        end
    else
        minimum(dataset[:, 2]) * 0.95, maximum(dataset[:, 2]) * 1.05
    end

    # Ensure Tmin < Tmax
    if Tmin > Tmax
        Tmin, Tmax = Tmax, Tmin
    end

    T_grid = temperature_grid(Tmin, Tmax, n_points)
    x = Vector{Float64}(undef, n_points)
    y = Vector{Float64}(undef, n_points)

    for (i, T) in enumerate(T_grid)
        state = interpolate_attractor_state(dataset, attractor_universe, τ, T)
        x[i] = xfn(state, cloud_slice)
        y[i] = yfn(state, cloud_slice)
    end

    return (; x, y)
end


