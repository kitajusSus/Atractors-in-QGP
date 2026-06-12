using LinearAlgebra
using NearestNeighbors
using Statistics

function dims(X; k = 12, tol = 1.0e-8)
    n, d = size(X)
    k = min(k, n)
    tree = KDTree(X')
    out = zeros(Float64, n)

    for i in 1:n
        idxs, _ = knn(tree, X[i, :], k, true)
        neighbors = X[idxs, :]

        μ = mean(neighbors, dims = 1)
        Y = neighbors .- μ

        C = (Y' * Y) / k
        λ = sort(eigvals(Symmetric(C)), rev = true)
        λ_norm = λ / sum(λ)

        out[i] = count(λ_norm .> tol)
    end

    return out
end


function normalizuj_max(X_features::AbstractMatrix{<:Real})
    max_vals = maximum(abs.(X_features), dims = 1)
    max_vals[max_vals .== 0.0] .= 1.0
    return X_features ./ max_vals
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
function compute_lpca(
        dataset_loadet,
        k_bazowe::Int,
        n_slices::Int;
        feature_cols::AbstractVector{<:Integer} = collect(2:size(dataset_loadet, 2))
    )
    taus = sort(unique(dataset_loadet[:, 1]))
    step = max(1, length(taus) ÷ n_slices)
    wybrane_tau = taus[1:step:end][1:n_slices]

    # k_zmienione = round(Int, k_bazowe * 2)
    tolerancja = 0.01

    # println("Tau | Śr. wymiar dla K=$(k_bazowe) | Śr. wymiar dla K=$(k_zmienione) | Zmiana (%)")

    # procenty = Float64[]
    tau_values = Float64[]
    mean_k1_values = Float64[]
    # mean_k2_values = Float64[]

    for tau in wybrane_tau
        idx, X_tau = get_tau_slice(dataset_loadet, tau; feature_cols = feature_cols)
        X_norm = normalizuj_max(X_tau)

        dims_k1 = dims(X_norm; k = k_bazowe, tol = tolerancja)
        # dims_k2 = dims(X_norm; k = k_zmienione, tol = tolerancja)

        mean_k1 = mean(dims_k1)
        # mean_k2 = mean(dims_k2)

        # diff_percent = mean_k1 > 0 ? ((mean_k2 - mean_k1) / mean_k1) * 100 : 0.0

        # println(tau, " | ", mean_k1, " | ", mean_k2, " | ", diff_percent, "%")
        # push!(procenty, diff_percent)
        push!(tau_values, tau)
        push!(mean_k1_values, mean_k1)
        # push!(mean_k2_values, mean_k2)
    end
    return tau_values, mean_k1_values #, mean_k2_values, procenty
end
