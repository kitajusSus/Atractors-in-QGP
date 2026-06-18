using LinearAlgebra
using NearestNeighbors
using Statistics

"""
    dims(X; k = 12, tol = 1e-8)

Lokalny wymiar danych wyznaczany metodą PCA.
"""
@views function dims(X; k = 12, tol = 1.0e-8)

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


"""
Normalizacja każdej cechy przez jej maksymalną wartość bezwzględną.
"""
function normalize_max(X)
    scale = maximum(abs.(X), dims = 1)
    for j in eachindex(scale)
        if scale[j] == 0
            scale[j] = 1
        end
    end
    return X ./ scale
end


"""
Generuje klasyczny swiss-roll.
"""
function swiss_roll(n; noise = 0.0)
    t = (3π / 2) .* (1 .+ 2 .* rand(n))
    h = 10 .* rand(n)
    x = t .* cos.(t)
    y = h
    z = t .* sin.(t)
    X = hcat(x, y, z)
    if noise > 0
        X .+= noise .* randn(size(X))
    end
    return X, t
end


"""
    compute_lpca(dataset, k, n_slices)

Średni lokalny wymiar danych jako funkcja czasu τ.
Zakładany format danych:

    [τ, x₁, x₂, ...]
"""
@views function compute_lpca(
        dataset,
        k::Int,
        n_slices::Int;
        feature_cols::AbstractVector{<:Integer} =
            collect(2:size(dataset, 2))
    )
    taus = sort(unique(dataset[:, 1]))

    idxs = round.(
        Int,
        range(1, length(taus), length = n_slices)
    )
    selected_taus = taus[idxs]

    tol = 0.01

    tau_values = Float64[]
    mean_dims = Float64[]
    std_dims = Float64[]

    for τ in selected_taus

        _, X_tau = get_tau_slice(
            dataset,
            τ;
            feature_cols = feature_cols
        )

        X_norm = normalize_max(X_tau)
        local_dims = dims(
            X_norm;
            k = k,
            tol = tol
        )

        push!(tau_values, τ)
        push!(mean_dims, mean(local_dims))
        push!(std_dims, std(local_dims))

    end

    return tau_values, mean_dims, std_dims
end


"""
    dynamic_lpca_analysis(...)

Analiza zapadania wymiaru lokalnego.
"""
@views function dynamic_lpca_analysis(
        dataset::AbstractMatrix{<:Real};
        K::Int = 15,
        eta::Real = 0.95,
        delta::Real = 0.05
    )

    taus = sort(unique(Float64.(dataset[:, 1])))
    d_bar = zeros(length(taus))
    feature_cols = 2:size(dataset, 2)
    for (nτ, τ) in enumerate(taus)

        mask = dataset[:, 1] .== τ
        X_tau = dataset[mask, feature_cols]
        n_points = size(X_tau, 1)
        k = min(K, n_points - 1)
        if k < 2
            d_bar[nτ] = 1.0
            continue
        end

        X_norm = normalize_max(X_tau)
        tree = KDTree(X_norm')

        idxs, _ = knn(
            tree,
            X_norm',
            k + 1,
            true
        )
        d_sum = 0
        for i in 1:n_points
            neighbors = idxs[i][2:end]
            X_local = X_norm[neighbors, :]
            X_local .-= mean(X_local, dims = 1)
            S = svdvals(X_local)
            λ = S .^ 2
            λ_sum = sum(λ)

            if λ_sum ≈ 0.0
                d_sum += 1
                continue
            end

            evr = 0.0
            d = 1
            for j in eachindex(λ)
                evr += λ[j] / λ_sum

                if evr >= eta
                    d = j
                    break
                end

            end
            d_sum += d
        end
        d_bar[nτ] = d_sum / n_points
    end
    tau_LPCA = NaN
    threshold = 1.0 + delta
    for i in eachindex(taus)
        if all(d_bar[i:end] .<= threshold)

            tau_LPCA = taus[i]

            break
        end
    end

    return (; taus, d_bar, tau_LPCA)
end


"""
    compute_lpca_entropy(...)

Entropia lokalnego PCA jako funkcja czasu τ.
"""
@views function compute_lpca_entropy(
        dataset::AbstractMatrix{<:Real},
        k::Int;
        n_slices::Union{Nothing, Int} = nothing,
        feature_cols::AbstractVector{<:Integer} =
            collect(2:size(dataset, 2)),
        tol::Real = 1.0e-8
    )

    taus = sort(unique(dataset[:, 1]))

    if n_slices === nothing
        selected_taus = taus
    else
        idxs = round.(
            Int,
            range(1, length(taus), length = n_slices)
        )

        selected_taus = taus[idxs]
    end

    # definicje tablic na dane
    tau_values = Float64[]
    mean_entropy = Float64[]
    std_entropy = Float64[]
    for τ in selected_taus
        _, X_tau = get_tau_slice(
            dataset,
            τ;
            feature_cols = feature_cols
        )
        X_norm = normalize_max(X_tau)
        n_points = size(X_norm, 1)
        k_local = min(k, n_points)

        if k_local < 2

            push!(tau_values, τ)
            push!(mean_entropy, 0.0)
            push!(std_entropy, 0.0)

            continue
        end
        tree = KDTree(X_norm')

        idxs, _ = knn(
            tree,
            X_norm',
            k_local,
            true
        )
        entropies = zeros(n_points)

        for i in 1:n_points
            neighbors = X_norm[idxs[i], :]
            μ = mean(neighbors, dims = 1)
            Y = neighbors .- μ
            C = (Y' * Y) / k_local
            λ = sort(
                eigvals(Symmetric(C)),
                rev = true
            )
            λ = max.(λ, 0.0)
            λ_sum = sum(λ)
            if λ_sum ≈ 0.0
                continue
            end
            λ_norm = λ / λ_sum
            r = count(λ_norm .> tol)

            if r <= 1
                continue
            end

            S = 0.0

            for val in λ_norm
                if val > tol
                    S -= val * log(val)
                end
            end
            entropies[i] = S / log(r)
        end
        push!(tau_values, τ)
        push!(mean_entropy, mean(entropies))
        push!(std_entropy, std(entropies))
    end
    return tau_values, mean_entropy, std_entropy
end


"""
    compute_stable_lpca_collapse(
        dataset::AbstractMatrix{<:Real};
        zakres_K::AbstractVector{<:Integer} = [10, 20, 30, 40],
        Scrit::Real = 0.2,
        delta_k::Real = 0.05,
        n_slices::Union{Nothing, Int} = nothing,
        feature_cols::AbstractVector{<:Integer} = collect(2:size(dataset, 2)),
        tol::Real = 1.0e-8,
        norm_type::Symbol = :embedding
    )

1. Mean spectral entropy over scales < Scrit
2. Standard deviation (spread) over scales < delta_k
"""
function compute_stable_lpca_collapse(
        dataset::AbstractMatrix{<:Real};
        zakres_K::AbstractVector{<:Integer} = [10, 20, 30, 40],
        Scrit::Real = 0.2,
        delta_k::Real = 0.05,
        n_slices::Union{Nothing, Int} = nothing,
        feature_cols::AbstractVector{<:Integer} = collect(2:size(dataset, 2)),
        tol::Real = 1.0e-8,
        norm_type::Symbol = :embedding
    )
    # Get unique proper times
    taus = sort(unique(dataset[:, 1]))
    if n_slices !== nothing
        idxs = round.(Int, range(1, length(taus), length = n_slices))
        wybrane_tau = taus[idxs]
    else
        wybrane_tau = taus
    end

    n_taus = length(wybrane_tau)
    n_k = length(zakres_K)
    k_max = maximum(zakres_K)

    S_PCA_matrix = zeros(Float64, n_taus, n_k)

    for (tau_idx, tau) in enumerate(wybrane_tau)
        _, X_tau = get_tau_slice(dataset, tau; feature_cols = feature_cols)
        X_norm = normalize_max(X_tau)

        n_points, d = size(X_norm)
        k_max_actual = min(k_max, n_points)
        if k_max_actual < 2
            continue
        end

        tree = KDTree(X_norm')
        knn_idxs, _ = knn(tree, X_norm', k_max_actual, true)

        for (k_idx, k_val) in enumerate(zakres_K)
            k_actual = min(k_val, n_points)
            if k_actual < 2
                S_PCA_matrix[tau_idx, k_idx] = 0.0
                continue
            end

            entropies = zeros(Float64, n_points)
            for i in 1:n_points
                neighbor_idxs = knn_idxs[i][1:k_actual]
                neighbors = X_norm[neighbor_idxs, :]

                μ = mean(neighbors, dims = 1)
                Y = neighbors .- μ

                C = (Y' * Y) / k_actual
                λ = sort(eigvals(Symmetric(C)), rev = true)
                λ = max.(λ, 0.0)
                sum_λ = sum(λ)

                if sum_λ ≈ 0.0
                    entropies[i] = 0.0
                    continue
                end

                λ_norm = λ / sum_λ
                S_i = 0.0
                for val in λ_norm
                    if val > tol
                        S_i -= val * log(val)
                    end
                end

                if norm_type == :embedding
                    entropies[i] = S_i / log(d)
                elseif norm_type == :active
                    r = count(λ_norm .> tol)
                    entropies[i] = r <= 1 ? 0.0 : S_i / log(r)
                else # :none
                    entropies[i] = S_i
                end
            end
            S_PCA_matrix[tau_idx, k_idx] = mean(entropies)
        end
    end

    S_PCA_mean = mean(S_PCA_matrix, dims = 2)[:]
    S_PCA_std = std(S_PCA_matrix, dims = 2)[:]

    tau_LPCA_k = fill(NaN, n_k)
    for (k_idx, k_val) in enumerate(zakres_K)
        for (tau_idx, tau) in enumerate(wybrane_tau)
            if S_PCA_matrix[tau_idx, k_idx] < Scrit
                tau_LPCA_k[k_idx] = tau
                break
            end
        end
    end

    tau_LPCA_stable = NaN
    for (tau_idx, tau) in enumerate(wybrane_tau)
        if S_PCA_mean[tau_idx] < Scrit && S_PCA_std[tau_idx] < delta_k
            tau_LPCA_stable = tau
            break
        end
    end

    return (;
        taus = wybrane_tau,
        S_PCA_matrix,
        S_PCA_mean,
        S_PCA_std,
        tau_LPCA_k,
        tau_LPCA_stable,
    )
end


"""
    compute_lpca_principal_angles(
        dataset::AbstractMatrix{<:Real};
        k::Int = 20,
        subspace_dim::Int = 1,
        n_slices::Union{Nothing, Int} = nothing,
        feature_cols::AbstractVector{<:Integer} = collect(2:size(dataset, 2))
    )

Compares neighboring local tangent spaces by computing the principal angles between them.
For each point x_i in a time slice, we find its nearest neighbor x_j (excluding itself).
We retrieve their local bases V_i and V_j (composed of the first `subspace_dim` principal components).
The cosines of the principal angles are the singular values of V_i' * V_j.
We return the angles in degrees: theta_r = arccos(sigma_r) * (180/pi).
"""
@views function compute_lpca_principal_angles(
        dataset::AbstractMatrix{<:Real};
        k::Int = 20,
        subspace_dim::Int = 1,
        n_slices::Union{Nothing, Int} = nothing,
        feature_cols::AbstractVector{<:Integer} = collect(2:size(dataset, 2))
    )
    taus = sort(unique(dataset[:, 1]))
    if n_slices !== nothing
        idxs = round.(Int, range(1, length(taus), length = n_slices))
        wybrane_tau = taus[idxs]
    else
        wybrane_tau = taus
    end

    mean_angles = Float64[]
    std_angles = Float64[]
    all_angles_per_tau = Vector{Vector{Float64}}()

    for tau in wybrane_tau
        _, X_tau = get_tau_slice(dataset, tau; feature_cols = feature_cols)
        X_norm = normalize_max(X_tau)

        n_points, d = size(X_norm)
        k_actual = min(k, n_points)
        if k_actual < 2
            push!(mean_angles, 0.0)
            push!(std_angles, 0.0)
            push!(all_angles_per_tau, Float64[])
            continue
        end

        tree = KDTree(X_norm')
        knn_idxs, _ = knn(tree, X_norm', k_actual, true)

        m = min(subspace_dim, d)
        V = Vector{Matrix{Float64}}(undef, n_points)

        for i in 1:n_points
            neighbor_idxs = knn_idxs[i]
            neighbors = X_norm[neighbor_idxs, :]
            μ = mean(neighbors, dims = 1)
            Y = neighbors .- μ

            F = svd(Y)
            V[i] = Matrix{Float64}(F.V[:, 1:m])
        end

        angles = Float64[]
        for i in 1:n_points
            if length(knn_idxs[i]) < 2
                continue
            end
            j = knn_idxs[i][2]

            M_proj = V[i]' * V[j]
            sigmas = svdvals(M_proj)
            sigmas = clamp.(sigmas, 0.0, 1.0)

            for σ in sigmas
                push!(angles, acos(σ) * (180.0 / π))
            end
        end

        if isempty(angles)
            push!(mean_angles, 0.0)
            push!(std_angles, 0.0)
            push!(all_angles_per_tau, Float64[])
        else
            push!(mean_angles, mean(angles))
            push!(std_angles, std(angles))
            push!(all_angles_per_tau, angles)
        end
    end

    return (;
        taus = wybrane_tau,
        mean_angles,
        std_angles,
        all_angles_per_tau,
    )
end
