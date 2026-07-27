using LinearAlgebra
using NearestNeighbors
using Statistics

function dims(X; k = 12, tol = 1.0e-8)
    n, d = size(X)
    k = min(k, n)
    tree = KDTree(X')
    macierz_zer = zeros(Float64, n)

    for i in 1:n
        idxs, _ = knn(tree, X[i, 2:4], k, true)
        neighbors = X[idxs, :]

        μ = mean(neighbors, dims = 1)
        Y = neighbors .- μ

        C = (Y' * Y) / k
        λ = sort(eigvals(Symmetric(C)), rev = true)
        λ_norm = λ / sum(λ)

        macierz_zer[i] = count(λ_norm .> tol)
    end
    return macierz_zer
end

function swiss_roll(n; noise = 0.0)
    t = (3π / 2) .* (1 .+ 2 .* rand(n))   # "angle"
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


# X, t = swiss_roll(2000, noise=0.1);
#
# dd = dims(X);
#
# meandim = mean(dd)
#
# print(meandim)
