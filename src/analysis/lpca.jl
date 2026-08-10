using LinearAlgebra
using NearestNeighbors
using Statistics


@views function dims(X; k = 12, tol = 1.0e-8)
    n, d = size(X)
    k = min(k, n)

    tree = KDTree(X')
    out = zeros(Float64, n)

    Y_buf = zeros(Float64, k, d)
    C_buf = zeros(Float64, d, d)

    knn_idxs, _ = knn(tree, X', k, true)

    for i in 1:n
        idxs = knn_idxs[i]

        Y = Y_buf[1:k, :]
        for r in 1:k
            for c in 1:d
                Y[r, c] = X[idxs[r], c]
            end
        end

        μ = mean(Y, dims = 1)
        Y .-= μ

        C = C_buf[1:d, 1:d]
        mul!(C, Y', Y)
        C ./= k

        λ = eigvals!(Symmetric(C))
        sum_λ = sum(λ)
        if sum_λ > 0.0
            cnt = 0
            for val in λ
                if (val / sum_λ) > tol
                    cnt += 1
                end
            end
            out[i] = cnt
        else
            out[i] = 1.0
        end
    end
    return out
end

@views function dims(X; k = 12, tol = 1.0e-8)
    n, d = size(X)
    k = min(k, n)

    tree = KDTree(X')
    out = zeros(Float64, n)

    Y_buf = zeros(Float64, k, d)
    C_buf = zeros(Float64, d, d)

    knn_idxs, _ = knn(tree, X', k, true)

    for i in 1:n
        idxs = knn_idxs[i]

        Y = Y_buf[1:k, :]
        for r in 1:k
            for c in 1:d
                Y[r, c] = X[idxs[r], c]
            end
        end

        μ = mean(Y, dims = 1)
        Y .-= μ

        C = C_buf[1:d, 1:d]
        mul!(C, Y', Y)
        C ./= k

        λ = eigvals!(Symmetric(C))
        sum_λ = sum(λ)
        if sum_λ > 0.0
            cnt = 0
            for val in λ
                if (val / sum_λ) > tol
                    cnt += 1
                end
            end
            out[i] = cnt
        else
            out[i] = 1.0
        end
    end
    return out
end

@views function dims(X; k = 12, tol = 1.0e-8)
    n, d = size(X)
    k = min(k, n)

    tree = KDTree(X')
    out = zeros(Float64, n)

    Y_buf = zeros(Float64, k, d)
    C_buf = zeros(Float64, d, d)

    knn_idxs, _ = knn(tree, X', k, true)

    for i in 1:n
        idxs = knn_idxs[i]

        Y = Y_buf[1:k, :]
        for r in 1:k
            for c in 1:d
                Y[r, c] = X[idxs[r], c]
            end
        end

        μ = mean(Y, dims = 1)
        Y .-= μ

        C = C_buf[1:d, 1:d]
        mul!(C, Y', Y)
        C ./= k

        λ = eigvals!(Symmetric(C))
        sum_λ = sum(λ)
        if sum_λ > 0.0
            cnt = 0
            for val in λ
                if (val / sum_λ) > tol
                    cnt += 1
                end
            end
            out[i] = cnt
        else
            out[i] = 1.0
        end
    end
    return out
end

function normalize_max(X)
    max_per_column = maximum(abs, X, dims = 1)
    max_per_column[max_per_column .== 0.0] .= 1.0
    return X ./ max_per_column
end

"""
    apply_normalization(X::AbstractMatrix{<:Real}, method::Union{Symbol, Function} = :max)

Applies data normalization / regularization to matrix `X`.

Supported methods:
- `:max`: Divides each column by max absolute value (rescales to [-1, 1]).
- `:minmax`: Min-Max column normalization (rescales to [0, 1]).
- `:zscore` / `:standard`: Zero mean, unit variance.
- `:none` / `:raw`: Unscaled raw matrix.
- Custom `Function`: User-supplied scaling function `f(X) -> X_scaled`.
"""
function apply_normalization(X::AbstractMatrix{<:Real}, method::Union{Symbol, Function} = :standard)
    if method isa Function
        return method(X)
    elseif method === :max
        return normalize_max(X)
    elseif method === :minmax
        Xn, _, _ = normalize_minmax(X)
        return Xn
    elseif method === :zscore || method === :standard
        Xf = Matrix{Float64}(X)
        μ = mean(Xf, dims = 1)
        σ = std(Xf, dims = 1)
        σ[σ .== 0.0] .= 1.0
        return (Xf .- μ) ./ σ
    elseif method === :none || method === :raw
        return Matrix{Float64}(X)
    else
        throw(ArgumentError("Unknown normalization method: $method. Supported: :max, :minmax, :zscore, :none, or a Function."))
    end
end

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

@views function compute_lpca(
        dataset,
        k::Int,
        n_slices::Int;
        feature_cols::AbstractVector{<:Integer} = collect(2:size(dataset, 2)),
        normalize::Union{Symbol, Function} = :standard
    )
    taus = sort(unique(dataset[:, 1]))
    idxs = round.(Int, range(1, length(taus), length = n_slices))
    selected_taus = taus[idxs]

    tol = 0.01
    n_selected = length(selected_taus)

    tau_values = zeros(Float64, n_selected)
    mean_dims = zeros(Float64, n_selected)
    std_dims = zeros(Float64, n_selected)

    for (idx, τ) in enumerate(selected_taus)
        _, X_tau = get_tau_slice(dataset, τ; feature_cols = feature_cols)
        X_norm = apply_normalization(X_tau, normalize)

        local_dims = dims(X_norm; k = k, tol = tol)

        tau_values[idx] = τ
        mean_dims[idx] = mean(local_dims)
        std_dims[idx] = std(local_dims)
    end

    return tau_values, mean_dims, std_dims
end


# nizej jest testowe nie ważnbe
@views function dynamic_lpca_analysis(
        dataset::AbstractMatrix{<:Real};
        K::Int = 15,
        eta::Real = 0.95,
        delta::Real = 0.05,
        feature_cols::Union{AbstractVector{<:Integer}, Nothing} = nothing
    )
    taus = sort(unique(Float64.(dataset[:, 1])))
    d_bar = zeros(length(taus))
    cols = isnothing(feature_cols) ? (2:size(dataset, 2)) : feature_cols
    d = length(cols)

    for (nτ, τ) in enumerate(taus)
        mask = dataset[:, 1] .== τ
        X_tau = dataset[mask, cols]
        n_points = size(X_tau, 1)
        k = min(K, n_points - 1)

        if k < 2
            d_bar[nτ] = 1.0
            continue
        end

        X_norm = normalize_max(X_tau)
        tree = KDTree(X_norm')

        idxs, _ = knn(tree, X_norm', k + 1, true)
        d_sum = 0

        X_local_buf = zeros(Float64, k, d)

        for i in 1:n_points
            neighbors = idxs[i][2:end]

            X_local = X_local_buf[1:k, :]
            for r in 1:k
                for c in 1:d
                    X_local[r, c] = X_norm[neighbors[r], c]
                end
            end

            X_local .-= mean(X_local, dims = 1)
            S = svdvals!(X_local)

            λ_sum = 0.0
            for s in S
                λ_sum += s^2
            end

            if λ_sum ≈ 0.0
                d_sum += 1
                continue
            end

            evr = 0.0
            d_val = 1
            for j in eachindex(S)
                evr += (S[j]^2) / λ_sum
                if evr >= eta
                    d_val = j
                    break
                end
            end
            d_sum += d_val
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


@views function compute_lpca_entropy(
        dataset::AbstractMatrix{<:Real},
        k::Int;
        n_slices::Union{Nothing, Int} = nothing,
        feature_cols::AbstractVector{<:Integer} = collect(2:size(dataset, 2)),
        tol::Real = 1.0e-8
    )
    taus = sort(unique(dataset[:, 1]))
    if n_slices === nothing
        selected_taus = taus
    else
        idxs = round.(Int, range(1, length(taus), length = n_slices))
        selected_taus = taus[idxs]
    end

    n_selected = length(selected_taus)
    tau_values = zeros(Float64, n_selected)
    mean_entropy = zeros(Float64, n_selected)
    std_entropy = zeros(Float64, n_selected)

    d = length(feature_cols)

    for (idx, τ) in enumerate(selected_taus)
        _, X_tau = get_tau_slice(dataset, τ; feature_cols = feature_cols)
        X_norm = normalize_max(X_tau)
        n_points = size(X_norm, 1)
        k_local = min(k, n_points)

        if k_local < 2
            tau_values[idx] = τ
            mean_entropy[idx] = 0.0
            std_entropy[idx] = 0.0
            continue
        end

        tree = KDTree(X_norm')
        idxs_knn, _ = knn(tree, X_norm', k_local, true)
        entropies = zeros(n_points)

        Y_buf = zeros(Float64, k_local, d)
        C_buf = zeros(Float64, d, d)

        for i in 1:n_points
            neighbors_idx = idxs_knn[i]

            Y = Y_buf[1:k_local, :]
            for r in 1:k_local
                for c in 1:d
                    Y[r, c] = X_norm[neighbors_idx[r], c]
                end
            end

            μ = mean(Y, dims = 1)
            Y .-= μ

            C = C_buf[1:d, 1:d]
            mul!(C, Y', Y)
            C ./= k_local

            λ = eigvals!(Symmetric(C))
            λ .= max.(λ, 0.0)
            λ_sum = sum(λ)

            if λ_sum ≈ 0.0
                continue
            end

            S = 0.0
            r_count = 0
            for val in λ
                λ_norm_val = val / λ_sum
                if λ_norm_val > tol
                    S -= λ_norm_val * log(λ_norm_val)
                    r_count += 1
                end
            end

            if r_count <= 1
                continue
            end

            entropies[i] = S / log(r_count)
        end

        tau_values[idx] = τ
        mean_entropy[idx] = mean(entropies)
        std_entropy[idx] = std(entropies)
    end
    return tau_values, mean_entropy, std_entropy
end

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
    d = length(feature_cols)

    S_PCA_matrix = zeros(Float64, n_taus, n_k)

    for (tau_idx, tau) in enumerate(wybrane_tau)
        _, X_tau = get_tau_slice(dataset, tau; feature_cols = feature_cols)
        X_norm = normalize_max(X_tau)

        n_points, d_size = size(X_norm)
        k_max_actual = min(k_max, n_points)
        if k_max_actual < 2
            continue
        end

        tree = KDTree(X_norm')
        knn_idxs, _ = knn(tree, X_norm', k_max_actual, true)

        Y_buf = zeros(Float64, k_max_actual, d)
        C_buf = zeros(Float64, d, d)

        for (k_idx, k_val) in enumerate(zakres_K)
            k_actual = min(k_val, n_points)
            if k_actual < 2
                S_PCA_matrix[tau_idx, k_idx] = 0.0
                continue
            end

            entropies = zeros(Float64, n_points)
            for i in 1:n_points
                @views neighbor_idxs = knn_idxs[i][1:k_actual]

                Y = Y_buf[1:k_actual, :]
                for r in 1:k_actual
                    for c in 1:d
                        Y[r, c] = X_norm[neighbor_idxs[r], c]
                    end
                end

                μ = mean(Y, dims = 1)
                Y .-= μ

                C = C_buf[1:d, 1:d]
                mul!(C, Y', Y)
                C ./= k_actual

                λ = eigvals!(Symmetric(C))
                λ .= max.(λ, 0.0)
                sum_λ = sum(λ)

                if sum_λ ≈ 0.0
                    entropies[i] = 0.0
                    continue
                end

                S_i = 0.0
                r_count = 0
                for val in λ
                    λ_norm = val / sum_λ
                    if λ_norm > tol
                        S_i -= λ_norm * log(λ_norm)
                        r_count += 1
                    end
                end

                if norm_type == :embedding
                    entropies[i] = S_i / log(d)
                elseif norm_type == :active
                    entropies[i] = r_count <= 1 ? 0.0 : S_i / log(r_count)
                else
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
