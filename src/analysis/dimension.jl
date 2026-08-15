using Statistics
using LinearAlgebra

"""
    estimate_effective_dimension(data::AbstractMatrix{<:Real})

Estymuje wymiar efektywny (participation ratio) z macierzy kowariancji dla DANYCH O DOWOLNEJ LICZBIE CECH.
"""
function estimate_effective_dimension(data::AbstractMatrix{<:Real})
    @assert size(data, 1) > 1 "Need at least two samples."
    X = Matrix{Float64}(data)
    Xc = X .- mean(X, dims = 1)
    C = cov(Xc)
    vals = eigvals(Symmetric(C))
    vals = real.(vals)
    vals = vals[vals .> 0]
    @assert !isempty(vals) "Covariance has no positive eigenvalues."
    return (sum(vals)^2) / sum(vals .^ 2)
end

const estimate_dimension = estimate_effective_dimension


"""
    scan_dimension_from_data(dataset::AbstractMatrix{<:Real}; atol=1.0e-3, feature_cols=nothing) -> NamedTuple

Skanuje dane w czasie i oblicza wymiar efektywny.
Jeśli `feature_cols=nothing`, automatycznie bierze pod uwagę WSZYSTKIE kolumny od drugiej do ostatniej.
"""
function scan_dimension_from_data(
        dataset::AbstractMatrix{<:Real};
        atol::Real = 1.0e-3,
        feature_cols::Union{AbstractVector{<:Integer}, Nothing} = nothing,
    )
    cols = isnothing(feature_cols) ? collect(2:size(dataset, 2)) : feature_cols
    @assert !isempty(cols) "Dane muszą mieć przynajmniej jedną kolumnę z cechami."

    taus = sort(unique(Float64.(dataset[:, 1])))
    valid_taus = Float64[]
    dimensions = Float64[]

    for tau in taus
        _, Xtau = get_tau_slice(dataset, tau; atol = atol, feature_cols = cols)

        if size(Xtau, 1) > size(Xtau, 2)
            Xtau_normalized, _, _ = normalize_minmax(Xtau)
            d = estimate_effective_dimension(Xtau_normalized)
            push!(valid_taus, tau)
            push!(dimensions, d)
        end
    end

    return (taus = valid_taus, dimensions = dimensions)
end


function to_2d_local_dimension(dataset::AbstractArray{<:Real, 3})
    dataset_2d = reshape(permutedims(dataset, (2, 1, 3)), :, size(dataset, 3))
    valid_rows = .!isnan.(dataset_2d[:, 1])
    return dataset_2d[valid_rows, :]
end


function estimate_effective_dimension(dataset::AbstractArray{<:Real, 3}, args...; kwargs...)
    return estimate_effective_dimension(to_2d_local_dimension(dataset), args...; kwargs...)
end
function scan_dimension_from_data(dataset::AbstractArray{<:Real, 3}, args...; kwargs...)
    return scan_dimension_from_data(to_2d_local_dimension(dataset), args...; kwargs...)
end
