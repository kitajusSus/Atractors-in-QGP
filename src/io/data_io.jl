using CSV
using DataFrames
using HDF5
using Serialization


function attractor_data(model::AbstractHydroModel; tau_max = 5.0, krok = 0.01, temperature_unit::Symbol = :fm)
    tspan_attr = (0.2, tau_max)

    ic_attr = [4, 7]
    sol_attr = solve_hydro(model, ic_attr, tspan_attr; saveat = krok)
    attractor_matrix = build_dataset([sol_attr]; temperature_unit = temperature_unit)

    return attractor_matrix
end

function attractor_data(
        model::MISModel;
        tau_max = 5.0,
        krok = 0.01,
        temperature_unit::Symbol = :fm,
        w_min::Real = 1.0e-4
    )
    w_max = Float64(tau_max)
    @assert w_min > 0 "w_min must be positive for the MIS attractor profile."
    @assert w_min < w_max "w_min must be smaller than tau_max."
    @assert krok > 0 "krok must be positive."

    p = model.params
    A0 = 6 * sqrt(p.eta_over_s / p.tau_pi)

    rhs_A(u, _, w) = begin
        A = u[1]
        dA = 18 * (8 * p.eta_over_s - w * A - (2 / 9) * p.tau_pi * A^2) /
            (p.tau_pi * w * (A + 12))
        return SVector(dA)
    end

    saveat = collect(Float64(w_min):Float64(krok):w_max)
    if saveat[end] < w_max
        push!(saveat, w_max)
    end

    problem = ODEProblem(rhs_A, SVector(Float64(A0)), (Float64(w_min), w_max), nothing)
    sol = solve(problem, Rodas5(); abstol = 1.0e-9, reltol = 1.0e-9, saveat = saveat)

    profile = Matrix{Float64}(undef, length(sol.t), 3)
    profile[:, 1] .= sol.t
    profile[:, 2] .= 1.0
    profile[:, 3] .= getindex.(sol.u, 1)
    return profile
end

"""
Save dataset matrix to CSV.
"""
function save_dataset_csv(path::AbstractString, data::AbstractMatrix{<:Real})
    @assert size(data, 2) >= 2 "Dataset must have at least columns [tau, feature1]."
    n_cols = size(data, 2)
    col_names = Symbol[]
    push!(col_names, :tau)
    if n_cols >= 2
        push!(col_names, :T)
    end
    if n_cols >= 3
        push!(col_names, :A)
    end
    for i in 4:n_cols
        push!(col_names, Symbol("feature$(i - 1)"))
    end
    df = DataFrame(data, col_names)
    CSV.write(path, df)
    return path
end

"""
Load dataset matrix from CSV.
"""
function load_dataset_csv(path::AbstractString)
    df = CSV.read(path, DataFrame)
    cols = names(df)

    name_map = Dict(lowercase(strip(String(c))) => c for c in cols)
    tau_col = get(name_map, "tau", nothing)
    t_col = get(name_map, "t", nothing)
    a_col = get(name_map, "a", nothing)

    if tau_col === nothing
        @assert size(df, 2) >= 2 "CSV must contain at least two columns."
        return Matrix{Float64}(df)
    end

    selected = [tau_col]
    if t_col !== nothing && t_col in cols
        push!(selected, t_col)
    end
    if a_col !== nothing && a_col in cols
        push!(selected, a_col)
    end

    for c in cols
        if !(c in selected)
            push!(selected, c)
        end
    end

    reordered_df = select(df, selected)
    return Matrix{Float64}(reordered_df)
end

"""
Save dataset matrix to HDF5.
"""
function save_dataset_h5(path::AbstractString, data::AbstractMatrix{<:Real})
    @assert size(data, 2) >= 2 "Dataset must have at least columns [tau, feature1]."
    h5open(path, "w") do f
        f["dataset"] = Matrix{Float64}(data)
    end
    return path
end

"""
Load dataset matrix from HDF5.
"""
function load_dataset_h5(path::AbstractString)
    return h5open(path, "r") do f
        @assert haskey(f, "dataset") "HDF5 file must contain /dataset"
        data = read(f["dataset"])
        @assert size(data, 2) >= 2 "Dataset must have at least columns [tau, feature1]."
        return Matrix{Float64}(data)
    end
end

"""
Save dataset using native Julia serialization (.jls).
"""
function save_dataset_jls(path::AbstractString, data::AbstractMatrix{<:Real})
    @assert size(data, 2) >= 2 "Dataset must have at least columns [tau, feature1]."
    open(path, "w") do io
        serialize(io, Matrix{Float64}(data))
    end
    return path
end

"""
Load dataset from native Julia serialization (.jls).
"""
function load_dataset_jls(path::AbstractString)
    return open(path, "r") do io
        data = deserialize(io)
        @assert data isa AbstractMatrix "Serialized object must be a matrix."
        @assert size(data, 2) >= 2 "Dataset must have at least columns [tau, feature1]."
        return Matrix{Float64}(data)
    end
end

"""
Save dataset by extension: .csv, .h5/.hdf5, .jls
"""
function save_dataset(path::AbstractString, data::AbstractMatrix{<:Real})
    lower = lowercase(path)
    if endswith(lower, ".csv")
        return save_dataset_csv(path, data)
    elseif endswith(lower, ".h5") || endswith(lower, ".hdf5")
        return save_dataset_h5(path, data)
    elseif endswith(lower, ".jls")
        return save_dataset_jls(path, data)
    else
        error("Unsupported format. Use .csv, .h5/.hdf5, or .jls")
    end
end

"""
Load dataset by extension: .csv, .h5/.hdf5, .jls
"""
function load_dataset(path::AbstractString)
    lower = lowercase(path)
    if endswith(lower, ".csv")
        return load_dataset_csv(path)
    elseif endswith(lower, ".h5") || endswith(lower, ".hdf5")
        return load_dataset_h5(path)
    elseif endswith(lower, ".jls")
        return load_dataset_jls(path)
    else
        error("Unsupported format. Use .csv, .h5/.hdf5, or .jls")
    end
end


"""
    attractor_state(attractor, τ, T, ncols)

Construct a virtual system state lying on the attractor.
"""
function attractor_state(
        attractor::AbstractMatrix{<:Real},
        τ::Real,
        T::Real,
        ncols::Int
    )
    # Wykrywanie jednostek temperatury (MeV vs fm^-1) i zestrojenie zakresów omega
    T_is_mev = T > 50.0
    T_attractor_is_mev = (sum(attractor[:, 2]) / size(attractor, 1)) > 50.0

    T_fm = T_is_mev ? T * FM_PER_MEV : T
    omega_target = τ * T_fm

    T_attractor_fm = T_attractor_is_mev ? attractor[:, 2] .* FM_PER_MEV : attractor[:, 2]
    omega_universe = attractor[:, 1] .* T_attractor_fm

    # Sortowanie pod kątem ekstrapolacji
    p_sort = sortperm(omega_universe)
    omega_sorted = omega_universe[p_sort]
    attr_sorted = attractor[p_sort, :]

    omega_min = omega_sorted[1]
    omega_max = omega_sorted[end]

    state = Vector{Float64}(undef, ncols)
    state[1] = τ
    state[2] = T

    if omega_target < omega_min
        # Ekstrapolacja dla małych omega: stała wartość początkowa
        row = attr_sorted[1, :]
        state[3] = row[3]
        if ncols > 3
            n_copy = min(ncols, size(attractor, 2))
            state[4:n_copy] .= row[4:n_copy]
            if ncols > n_copy
                state[(n_copy + 1):end] .= 0.0
            end
        end
    elseif omega_target > omega_max
        # Ekstrapolacja dla dużych omega: spadek Navier-Stokesa A ~ 1/omega
        A_max = attr_sorted[end, 3]
        state[3] = A_max * (omega_max / omega_target)
        if ncols > 3
            row = attr_sorted[end, :]
            n_copy = min(ncols, size(attractor, 2))
            state[4:n_copy] .= row[4:n_copy]
            if ncols > n_copy
                state[(n_copy + 1):end] .= 0.0
            end
        end
    else
        # Dopasowanie wewnątrz zakresu
        idx = argmin(abs.(omega_sorted .- omega_target))
        row = attr_sorted[idx, :]
        state[3] = row[3]
        if ncols > 3
            n_copy = min(ncols, size(attractor, 2))
            state[4:n_copy] .= row[4:n_copy]
            if ncols > n_copy
                state[(n_copy + 1):end] .= 0.0
            end
        end
    end

    return state
end


"""
    get_attractor_line(
        dataset,
        attractor,
        τ,
        xdef,
        ydef;
        n_points=150
    )

Return coordinates of the attractor line projected onto
the chosen observables.
"""
function get_attractor_line(
        dataset::AbstractMatrix{<:Real},
        attractor::AbstractMatrix{<:Real},
        τ::Real,
        xdef,
        ydef;
        n_points::Int = 150
    )

    _, xmap = resolve_def(xdef)
    _, ymap = resolve_def(ydef)

    _, slice = get_tau_slice(dataset, τ; feature_cols = 1:size(dataset, 2))

    Tmin = minimum(slice[:, 2])
    Tmax = maximum(slice[:, 2])

    Tgrid = range(
        Tmin * 0.95,
        Tmax * 1.05;
        length = n_points
    )

    x = Vector{Float64}(undef, n_points)
    y = Vector{Float64}(undef, n_points)

    for (i, T) in enumerate(Tgrid)

        state = attractor_state(
            attractor,
            τ,
            T,
            size(dataset, 2)
        )

        x[i] = xmap(state, slice)
        y[i] = ymap(state, slice)

    end

    return (; x, y)
end
