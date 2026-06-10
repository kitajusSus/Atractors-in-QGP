"""
    jacobian_analysis.jl  —  Jacobian-based dimensionality-reduction analysis for PINN

Implements the method:

    "Dimensionality reduction via parametric PINN + automatic differentiation"

The PINN is treated as a **global**, differentiable model of phase space:

    x̂_θ(τ ; x₀) = (T̂(τ), Â(τ))

where x₀ = (T₀, A₀, η/s, τ_π, …) are initial / parameter coordinates.

At each proper time τ we compute the Jacobian

    J(τ) = ∂x̂_θ / ∂x₀   ∈ ℝ^{n_out × n_ic}

using ForwardDiff (forward-mode AD) — exact (machine-precision), no finite differences.

The singular values {σ_i} of J define the **effective (participation-ratio) dimension**

    d_eff(τ) = (Σ σ_i²)² / Σ σ_i⁴  ∈ [1, n_out]

Early times: d_eff ≈ n_out  (system remembers initial conditions)
Late times:  d_eff → 1      (attractor — one-dimensional manifold)

For higher-dimensional problems every quantity generalises naturally via SVD.

Public API
──────────
    pinn_jacobian(result, τ, T₀, A₀, η, τπ)           → Matrix{Float64}  (2×2)
    pinn_jacobian_full(result, τ, T₀, A₀, η, τπ)       → Matrix{Float64}  (2×4)
    d_eff_from_singular_values(σ)                        → Float64
    pinn_deff_at(result, τ, T₀, A₀, η, τπ; full)        → Float64
    pinn_deff_scan(result, taus, ic_ensemble; ...)        → NamedTuple
    pinn_jacobian_scan(result, taus, ic_ensemble; ...)    → NamedTuple
    pinn_hydrodynamisation_time(scan; threshold)          → Float64
    sample_ic_ensemble(result, n; seed)                   → Vector
    fixed_transport_ic_ensemble(result, n; η, τπ, seed)   → Vector
    pinn_dimensionality_workflow(model, config; ...)       → NamedTuple
"""

using ForwardDiff
using LinearAlgebra
using Statistics
using ProgressMeter

# ─────────────────────────────────────────────────────────────────────────────
# Internal helpers
# ─────────────────────────────────────────────────────────────────────────────

"""
    _pinn_forward(nn, ps, st, τ_n, ic_n, config)

Raw normalised forward pass.
ic_n : 4-element vector [T₀_n, A₀_n, η_n, τπ_n] in normalised space.
Returns the 2-element normalised output (tracked by ForwardDiff dual numbers).
"""
function _pinn_forward(nn, ps, st, τ_n::Real, ic_n::AbstractVector, config::PINNConfig)
    x = vcat(τ_n, ic_n)          # 5-dim normalised input
    y, _ = nn(x, ps, st)
    return y                       # 2-element normalised output
end

"""
    _build_ic_norm(T₀, A₀, η, τπ, config) → Vector{Float64}

Normalise the 4-element IC block (everything except τ) to [-1,1].
"""
@inline function _build_ic_norm(T₀::Real, A₀::Real, η::Real, τπ::Real, cfg::PINNConfig)
    return Float64[
        _norm(T₀, cfg.T_range...),
        _norm(A₀, cfg.A_range...),
        _norm(η, cfg.η_range...),
        _norm(τπ, cfg.τπ_range...),
    ]
end

# ─────────────────────────────────────────────────────────────────────────────
# Core: Jacobian computation
# ─────────────────────────────────────────────────────────────────────────────

"""
    pinn_jacobian(result::PINNResult, τ, T₀, A₀, η, τπ) → Matrix{Float64}  (2×2)

Exact Jacobian  J(τ) = ∂(T̂, Â)/∂(T₀, A₀)  via ForwardDiff.

Transport coefficients (η, τπ) are held fixed.
Result is in physical units (fm⁻¹ for T, dimensionless for A).
"""
function pinn_jacobian(
        result::PINNResult,
        τ::Real, T₀::Real, A₀::Real, η::Real, τπ::Real,
    )
    nn, ps, st, cfg = result.network, result.ps, result.st, result.config

    τ_n = _norm(τ, cfg.τ_range...)
    Tlo, Thi = cfg.T_range
    Alo, Ahi = cfg.A_range

    # Physical output scale  (d_phys / d_norm)
    scale_out = Float64[(Thi - Tlo) / 2, (Ahi - Alo) / 2]

    # Physical input scale   (d_norm / d_phys)
    inv_scale_in = Float64[2.0 / (Thi - Tlo), 2.0 / (Ahi - Alo)]

    ic_n_fixed = _build_ic_norm(T₀, A₀, η, τπ, cfg)

    # Differentiate only w.r.t. [T₀_norm, A₀_norm]
    f = x0n -> begin
        ic_full = vcat(x0n, ic_n_fixed[3:end])   # keep η_n, τπ_n fixed
        _pinn_forward(nn, ps, st, τ_n, ic_full, cfg)
    end

    J_norm = ForwardDiff.jacobian(f, ic_n_fixed[1:2])   # 2×2 normalised

    # Convert:  J_phys[i,j] = scale_out[i] * J_norm[i,j] * inv_scale_in[j]
    J_phys = J_norm .* (scale_out .* inv_scale_in')
    return J_phys
end


"""
    pinn_jacobian_full(result::PINNResult, τ, T₀, A₀, η, τπ) → Matrix{Float64}  (2×4)

Full 2×4 Jacobian:  ∂(T̂, Â) / ∂(T₀, A₀, η, τπ).

Useful for multidimensional sensitivity analysis where transport parameters
are also treated as uncertain / variable.
"""
function pinn_jacobian_full(
        result::PINNResult,
        τ::Real, T₀::Real, A₀::Real, η::Real, τπ::Real,
    )
    nn, ps, st, cfg = result.network, result.ps, result.st, result.config

    τ_n = _norm(τ, cfg.τ_range...)
    Tlo, Thi = cfg.T_range;  Alo, Ahi = cfg.A_range
    ηlo, ηhi = cfg.η_range;  τπlo, τπhi = cfg.τπ_range

    scale_out = Float64[(Thi - Tlo) / 2, (Ahi - Alo) / 2]
    inv_scale_in = Float64[2 / (Thi - Tlo), 2 / (Ahi - Alo), 2 / (ηhi - ηlo), 2 / (τπhi - τπlo)]

    f = ic_n -> _pinn_forward(nn, ps, st, τ_n, ic_n, cfg)

    x0_norm = _build_ic_norm(T₀, A₀, η, τπ, cfg)
    J_norm = ForwardDiff.jacobian(f, x0_norm)   # 2×4

    J_phys = J_norm .* (scale_out .* inv_scale_in')
    return J_phys
end

# ─────────────────────────────────────────────────────────────────────────────
# d_eff: participation-ratio dimension
# ─────────────────────────────────────────────────────────────────────────────

"""
    d_eff_from_singular_values(σ) → Float64

Effective (participation-ratio) dimension:

    d_eff = (Σ σᵢ²)² / Σ σᵢ⁴

= n_out when all σ equal (isotropic phase space)
= 1     when only one σ non-zero (collapsed to attractor)

Returns NaN if all singular values are zero.
"""
function d_eff_from_singular_values(σ::AbstractVector{<:Real})
    λ = σ .^ 2
    s1 = sum(λ)
    s2 = sum(λ .^ 2)
    (s2 <= 0.0 || s1 <= 0.0) && return NaN
    return s1^2 / s2
end


"""
    pinn_deff_at(result, τ, T₀, A₀, η, τπ; full=false) → Float64

d_eff at a single (τ, x₀) point.
`full=true` → use 2×4 Jacobian (all four IC parameters).
"""
function pinn_deff_at(
        result::PINNResult,
        τ::Real, T₀::Real, A₀::Real, η::Real, τπ::Real;
        full::Bool = false,
    )
    J = full ?
        pinn_jacobian_full(result, τ, T₀, A₀, η, τπ) :
        pinn_jacobian(result, τ, T₀, A₀, η, τπ)
    return d_eff_from_singular_values(svdvals(J))
end

# ─────────────────────────────────────────────────────────────────────────────
# d_eff scan over ensemble
# ─────────────────────────────────────────────────────────────────────────────

"""
    pinn_deff_scan(result, taus, ic_ensemble; full, verbose) → NamedTuple

Scan d_eff(τ) over a time grid and an IC ensemble.
Each element of `ic_ensemble` must be indexable as ic[1..4] = (T₀, A₀, η, τπ).

Returns:
- `taus`         : time grid
- `d_eff_mean`   : ensemble-mean d_eff(τ)
- `d_eff_std`    : ensemble-std  d_eff(τ)
- `d_eff_matrix` : N_tau × N_ic raw values
- `sigma1_mean`  : mean σ₁(τ)  (dominant singular value)
- `sigma2_mean`  : mean σ₂(τ)  (suppressed singular value)
- `sigma_ratio`  : mean σ₂/σ₁  (→ 0 at hydrodynamisation)
"""
function pinn_deff_scan(
        result::PINNResult,
        taus::AbstractVector{<:Real},
        ic_ensemble;
        full::Bool = false,
        verbose::Bool = true,
    )
    N_tau = length(taus)
    N_ic = length(ic_ensemble)

    d_mat = fill(NaN, N_tau, N_ic)
    σ1_mat = fill(NaN, N_tau, N_ic)
    σ2_mat = fill(NaN, N_tau, N_ic)

    prog = verbose ?
        Progress(N_tau * N_ic; desc = "Jacobian d_eff scan: ", showspeed = true) :
        nothing

    for (j, ic) in enumerate(ic_ensemble)
        T₀ = Float64(ic[1]);  A₀ = Float64(ic[2])
        η = Float64(ic[3]);  τπ = Float64(ic[4])

        for (i, τ) in enumerate(taus)
            try
                J = full ?
                    pinn_jacobian_full(result, τ, T₀, A₀, η, τπ) :
                    pinn_jacobian(result, τ, T₀, A₀, η, τπ)
                σ = svdvals(J)
                d_mat[i, j] = d_eff_from_singular_values(σ)
                σ1_mat[i, j] = length(σ) >= 1 ? σ[1] : NaN
                σ2_mat[i, j] = length(σ) >= 2 ? σ[2] : NaN
            catch
                # leave NaN on failure (e.g. ODE domain error)
            end
            isnothing(prog) || next!(prog)
        end
    end

    _nm(M) = Float64[
        let v = filter(isfinite, M[i, :])
                isempty(v) ? NaN : mean(v)
        end
            for i in 1:N_tau
    ]
    _ns(M) = Float64[
        let v = filter(isfinite, M[i, :])
                length(v) > 1 ? std(v) : 0.0
        end
            for i in 1:N_tau
    ]

    return (
        taus = collect(Float64, taus),
        d_eff_mean = _nm(d_mat),
        d_eff_std = _ns(d_mat),
        d_eff_matrix = d_mat,
        sigma1_mean = _nm(σ1_mat),
        sigma2_mean = _nm(σ2_mat),
        sigma_ratio = _nm(σ2_mat ./ max.(σ1_mat, fill(1.0e-30, size(σ1_mat)))),
    )
end

# ─────────────────────────────────────────────────────────────────────────────
# Full Jacobian scan (returns raw J matrices)
# ─────────────────────────────────────────────────────────────────────────────

"""
    pinn_jacobian_scan(result, taus, ic_ensemble; full, verbose) → NamedTuple

Raw Jacobian matrices and singular values for every (τ, x₀) pair.

Returns:
- `taus`       : time grid
- `ics`        : ic ensemble
- `jacobians`  : N_tau × N_ic Matrix{Matrix{Float64}}
- `svd_values` : N_tau × N_ic Matrix{Vector{Float64}}
"""
function pinn_jacobian_scan(
        result::PINNResult,
        taus::AbstractVector{<:Real},
        ic_ensemble;
        full::Bool = false,
        verbose::Bool = true,
    )
    N_tau = length(taus)
    N_ic = length(ic_ensemble)

    Js = Matrix{Matrix{Float64}}(undef, N_tau, N_ic)
    svds = Matrix{Vector{Float64}}(undef, N_tau, N_ic)

    prog = verbose ?
        Progress(N_tau * N_ic; desc = "Jacobian full scan: ", showspeed = true) :
        nothing

    for (j, ic) in enumerate(ic_ensemble)
        T₀ = Float64(ic[1]);  A₀ = Float64(ic[2])
        η = Float64(ic[3]);  τπ = Float64(ic[4])
        for (i, τ) in enumerate(taus)
            J = full ?
                pinn_jacobian_full(result, τ, T₀, A₀, η, τπ) :
                pinn_jacobian(result, τ, T₀, A₀, η, τπ)
            Js[i, j] = J
            svds[i, j] = svdvals(J)
            isnothing(prog) || next!(prog)
        end
    end

    return (
        taus = collect(Float64, taus), ics = ic_ensemble,
        jacobians = Js, svd_values = svds,
    )
end


"""
    pinn_hydrodynamisation_time(scan; threshold=1.05) → Float64

First τ where d_eff_mean < threshold.  Returns NaN if never reached.
"""
function pinn_hydrodynamisation_time(scan::NamedTuple; threshold::Real = 1.05)
    idx = findfirst(d -> isfinite(d) && d < threshold, scan.d_eff_mean)
    isnothing(idx) && return NaN
    return scan.taus[idx]
end


"""
    sample_ic_ensemble(result::PINNResult, n; seed=7) → Vector

Sample n ICs (T₀, A₀, η, τπ) uniformly from `result.config` ranges.
"""
function sample_ic_ensemble(result::PINNResult, n::Int)
    rng = Xoshiro(5)
    cfg = result.config
    _r(lo, hi) = lo + rand(rng) * (hi - lo)
    return [
        (
                T₀ = _r(cfg.T_range...), A₀ = _r(cfg.A_range...),
                η = _r(cfg.η_range...), τπ = _r(cfg.τπ_range...),
            )
            for _ in 1:n
    ]
end

"""
    fixed_transport_ic_ensemble(result, n; η, τπ, seed=7) → Vector

Like `sample_ic_ensemble` but with fixed transport coefficients.
"""
function fixed_transport_ic_ensemble(
        result::PINNResult, n::Int;
        η::Real, τπ::Real
    )
    rng = Xoshiro(5)
    cfg = result.config
    _r(lo, hi) = lo + rand(rng) * (hi - lo)
    return [
        (
                T₀ = _r(cfg.T_range...), A₀ = _r(cfg.A_range...),
                η = Float64(η), τπ = Float64(τπ),
            )
            for _ in 1:n
    ]
end

# ─────────────────────────────────────────────────────────────────────────────
# Full pipeline
# ─────────────────────────────────────────────────────────────────────────────

"""
    pinn_dimensionality_workflow(model, config=PINNConfig(); kwargs...) → NamedTuple

End-to-end pipeline:
1. Train PINN
2. Sample IC ensemble
3. d_eff scan over taus
4. Detect hydrodynamisation time

Keyword arguments:
- `n_ic`       : ensemble size (default 200)
- `taus`       : time grid (default 0.22:0.02:1.2)
- `full`       : use full 2×4 Jacobian (default false)
- `seed_ic`    : seed for IC sampling (default 7)
- `train_kw…`  : forwarded to `train_pinn` (n_epochs, learning_rate, etc.)

Returns `(result, ics, scan, τ_hyd)`.
"""
function pinn_dimensionality_workflow(
        model::AbstractHydroModel,
        config::PINNConfig = PINNConfig();
        n_ic::Int = 200,
        taus = collect(0.22:0.02:1.2),
        full::Bool = false,
        train_kw...
    )
    result = train_pinn(model, config; train_kw...)
    ics = sample_ic_ensemble(result, n_ic)
    scan = pinn_deff_scan(result, taus, ics; full = full)
    τ_hyd = pinn_hydrodynamisation_time(scan)
    return (result = result, ics = ics, scan = scan, τ_hyd = τ_hyd)
end
