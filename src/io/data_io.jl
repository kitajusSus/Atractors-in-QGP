using CSV
using DataFrames
using HDF5
using Serialization

"""
Save dataset matrix to CSV.
"""
function save_dataset_csv(path::AbstractString, data::AbstractArray{<:Real})
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
function save_dataset_h5(path::AbstractString, data::AbstractArray{<:Real})
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
function save_dataset_jls(path::AbstractString, data::AbstractArray{<:Real})
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
function save_dataset(path::AbstractString, data::AbstractArray{<:Real})
    flat_data = data isa AbstractArray{<:Real, 3} ? reshape(data, :, size(data, 3)) : data
    lower = lowercase(path)
    if endswith(lower, ".csv")
        return save_dataset_csv(path, flat_data)
    elseif endswith(lower, ".h5") || endswith(lower, ".hdf5") || endswith(lower, ".hd5")
        return save_dataset_h5(path, flat_data)
    elseif endswith(lower, ".jls")
        return save_dataset_jls(path, flat_data)
    else
        error("Unsupported format. Use .csv, .h5/.hdf5/.hd5, or .jls")
    end
end

"""
Load dataset by extension: .csv, .h5/.hdf5, .jls
"""
function load_dataset(path::AbstractString)
    lower = lowercase(path)
    if endswith(lower, ".csv")
        return load_dataset_csv(path)
    elseif endswith(lower, ".h5") || endswith(lower, ".hdf5") || endswith(lower, ".hd5")
        return load_dataset_h5(path)
    elseif endswith(lower, ".jls")
        return load_dataset_jls(path)
    else
        error("Unsupported format. Use .csv, .h5/.hdf5/.hd5, or .jls")
    end
end
