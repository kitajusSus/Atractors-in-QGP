using LinearAlgebra
using NearestNeighbors
using Statistics

function dims(X; k=12, tol=0.1)
    # X of size (n_samples, n_features).
    n, d = size(X)

    # Build KD-tree for neighbor search
    # NearestNeighbors expects columns = points
    tree = KDTree(permutedims(X))  

    dims = zeros(Int, n)

    for i in 1:n
        # Find k nearest neighbors (including itself)
        idxs, _ = knn(tree, X[i, :], k, true)

        # Extract neighbors
        neighbors = X[idxs, :]

        # Center the neighborhood
        μ = mean(neighbors, dims=1)
        Y = neighbors .- μ

        # Covariance matrix
        C = (Y' * Y) / k

        # Eigenvalues (sorted descending)
        λ = sort(eigvals(C), rev=true)

        # Normalize eigenvalues
        λ_norm = λ / sum(λ)

        # Estimate dimension: count significant components
        dims[i] = count(λ_norm .> tol)
    end

    return dims
end

using Random

function swiss_roll(n; noise=0.0)
    t = (3π/2) .* (1 .+ 2 .* rand(n))   # "angle"
    h = 10 .* rand(n)                   # height

    x = t .* cos.(t)
    y = h
    z = t .* sin.(t)

    X = hcat(x, y, z)

    if noise > 0
        X .+= noise .* randn(size(X))
    end

    return X, t
end


X, t = swiss_roll(2000, noise=0.1);

dd = dims(X);

meandim = mean(dd)

print(meandim)

