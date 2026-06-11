using LinearAlgebra
using KrylovKit
using Statistics
using NearestNeighbors
using SparseArrays


"""
Locally Linear Embedding (LLE)
X: matrix of input data with dimensions (D, N), `D` dimension of features,  `N` number of samples
k: numer of nearest neighbors to consider for each point
d: target dimension of the embedding space
"""
function lle(X::AbstractMatrix{T}; k::Int = 10, d::Int = 2, kdtree = nothing) where {T <: AbstractFloat}
    D_dim, N = size(X)

    I_idx = Int[]
    J_idx = Int[]
    V = T[]
    sizehint!(I_idx, N * k)
    sizehint!(J_idx, N * k)
    sizehint!(V, N * k)

    if isnothing(kdtree)
        kdtree = KDTree(X)
    end

    idxs, _ = knn(kdtree, X, k + 1)

    for i in 1:N
        neighbors = @view idxs[i][2:(k + 1)]

        Z = X[:, neighbors] .- X[:, i]
        C = Z' * Z
        C += I(length(neighbors)) * 1.0e-3 * tr(C)

        w = C \ ones(T, length(neighbors))
        w ./= sum(w)

        append!(I_idx, fill(i, length(neighbors)))
        append!(J_idx, neighbors)
        append!(V, w)
    end

    W = sparse(I_idx, J_idx, V, N, N)

    M = (I - W)' * (I - W)

    sigma = -1.0e-5
    F = cholesky(Symmetric(M - sigma * I))
    inv_f = x -> F \ x

    _, vecs, _ = eigsolve(inv_f, N, d + 1, :LM; ishermitian = true, maxiter = 300, tol = 1.0e-5)

    return reduce(hcat, vecs[2:(d + 1)])
end

function run_lle_per_time(
        dataset::AbstractMatrix{<:Real};
        k::Integer = 10,
        d::Integer = 2,
        atol::Real = 1.0e-8,
        feature_cols::AbstractVector{<:Integer} = collect(2:size(dataset, 2))
    )

    taus = sort(unique(Float64.(dataset[:, 1])))
    lle_results = Dict{Float64, Matrix{Float64}}()

    for tau in taus
        _, Xtau = get_tau_slice(dataset, tau; atol = atol, feature_cols = feature_cols)

        if size(Xtau, 1) > k
            X_input = Matrix{Float64}(Xtau)'
            res = lle(X_input; k = k, d = d)
            lle_results[tau] = res
        end
    end

    return (taus = taus, lle_results = lle_results)
end


function run_lle_for_selected_taus(
        dataset::AbstractMatrix{<:Real},
        target_taus::Vector{Float64};
        k::Integer = 20,
        d::Integer = 2,
        atol::Real = 1.0e-3,
        feature_cols::AbstractVector{<:Integer} = collect(2:size(dataset, 2))
    )
    lle_results = Dict{Float64, Matrix{Float64}}()

    for tau in target_taus
        _, Xtau = get_tau_slice(dataset, tau; atol = atol, feature_cols = feature_cols)

        if size(Xtau, 1) > k
            println("Obliczam LLE dla tau = $tau (liczba punktów: $(size(Xtau, 1)))...")
            res = lle(Matrix{Float64}(Xtau)'; k = k, d = d)
            lle_results[tau] = Matrix(res')
        end
    end

    return (taus = target_taus, lle_results = lle_results)
end


function lle_spectrum(X; k = 10, kdtree = nothing, nλ = 10)
    D, N = size(X)

    I_idx = Int[]
    J_idx = Int[]
    V = Float64[]
    sizehint!(I_idx, N * k)
    sizehint!(J_idx, N * k)
    sizehint!(V, N * k)

    if isnothing(kdtree)
        kdtree = KDTree(X)
    end

    idxs, _ = knn(kdtree, X, k + 1)

    for i in 1:N
        neighbor_idx = @view idxs[i][2:(k + 1)]

        Z = X[:, neighbor_idx] .- X[:, i]
        C = Z' * Z + I(length(neighbor_idx)) * 1.0e-3 * tr(Z' * Z)

        w = C \ ones(length(neighbor_idx))
        w ./= sum(w)

        append!(I_idx, fill(i, length(neighbor_idx)))
        append!(J_idx, neighbor_idx)
        append!(V, w)
    end

    W = sparse(I_idx, J_idx, V, N, N)

    M = (I - W)' * (I - W)

    sigma = -1.0e-5
    F = cholesky(Symmetric(M - sigma * I))
    inv_f = x -> F \ x

    vals_inv, _, _ = eigsolve(inv_f, N, nλ, :LM; ishermitian = true, maxiter = 300, tol = 1.0e-5)

    vals = (1 ./ real(vals_inv)) .+ sigma
    return sort(vals)
end


function lle_spectrum_over_k(dataset; tau, k_values = 5:5:50, atol = 1.0e-8, nλ = 10)
    _, Xτ = get_tau_slice(dataset, tau; atol = atol)

    X = Matrix{Float64}(Xτ)'

    spectra = Vector{Any}(undef, length(k_values))
    tree = KDTree(X)

    Threads.@threads for i in 1:length(k_values)
        k = k_values[i]
        if size(X, 2) > k + 2
            spectra[i] = lle_spectrum(X; k = k, kdtree = tree, nλ = nλ)
        end
    end

    valid_mask = isassigned.(Ref(spectra), 1:length(spectra))

    return spectra[valid_mask], k_values[valid_mask]
end

"""
Liczy średnią i odchylenie standardowe dla spektrum wartości własnych 
lub czegokolwiek innego 
"""
function spectrum_statistics(spectra)
    max_len = minimum(length.(spectra))
    S = reduce(hcat, [s[1:max_len] for s in spectra])

    μ = mean(S, dims = 2)[:]
    σ = std(S, dims = 2)[:]

    return μ, σ
end


function scan_lle_spectrum(dataset; taus, k_values = 5:5:50)
    results = Dict()

    for τ in taus
        spectra, ks = lle_spectrum_over_k(dataset; tau = τ, k_values = k_values)

        if length(spectra) < 2
            continue
        end

        μ, σ = spectrum_statistics(spectra)

        results[τ] = (mean = μ, std = σ, k = ks)
    end

    return results
end
