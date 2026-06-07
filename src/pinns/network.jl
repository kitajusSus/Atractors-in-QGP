using Lux
using Random
using ComponentArrays

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

"""
    PINNConfig

Hyperparameters and normalization ranges for the PINN.

All physical variables are normalized to [-1, 1] before being fed to the
network. The ranges should bracket the full space of interest.

# Fields
- `hidden_layers`: number of hidden layers (not counting input/output)
- `hidden_size`:   neurons per hidden layer
- `activation`:    activation function — tanh recommended for ODEs
- `τ_range`:       (τ_min, τ_max) in fm/c
- `T_range`:       (T_min, T_max) in fm⁻¹  (internal unit)
- `A_range`:       (A_min, A_max) dimensionless
- `η_range`:       (η/s min, η/s max)
- `τπ_range`:      (τ_π min, τ_π max) in fm/c
"""
struct PINNConfig{F}
    hidden_layers::Int
    hidden_size::Int
    activation::F
    τ_range::Tuple{Float64, Float64}
    T_range::Tuple{Float64, Float64}
    A_range::Tuple{Float64, Float64}
    η_range::Tuple{Float64, Float64}
    τπ_range::Tuple{Float64, Float64}
end

"""
    PINNConfig(; kwargs...)

Create a `PINNConfig` with sensible defaults for Bjorken flow in fm units.

Default ranges cover the physically relevant region for BRSSS/MIS models.
"""
function PINNConfig(;
        hidden_layers::Int = 4,
        hidden_size::Int = 64,
        activation = tanh,
        τ_range = (0.1, 2.0),
        T_range = (0.5, 15.0),   # fm⁻¹  (~100–3000 MeV)
        A_range = (-15.0, 25.0),
        η_range = (0.05, 1.0),
        τπ_range = (0.05, 1.0),
    )
    F = typeof(activation)
    return PINNConfig{F}(
        hidden_layers, hidden_size, activation,
        Float64.(τ_range), Float64.(T_range),
        Float64.(A_range), Float64.(η_range),
        Float64.(τπ_range),
    )
end

# ---------------------------------------------------------------------------
# Network construction
# ---------------------------------------------------------------------------

"""
    build_pinn_network(config::PINNConfig) -> Lux.Chain

Build a fully-connected MLP for the PINN.

**Input** (5-dim, normalized to [-1,1]):
    [τ_n, T₀_n, A₀_n, (η/s)_n, τ_π_n]

**Output** (2-dim, normalized):
    [T_n(τ), A_n(τ)]

Denormalize outputs with `denormalize_pinn_output` to get physical values.
"""
function build_pinn_network(config::PINNConfig)
    hidden = ntuple(
        _ -> Dense(config.hidden_size, config.hidden_size, config.activation),
        config.hidden_layers - 1,
    )
    return Chain(
        Dense(5, config.hidden_size, config.activation),
        hidden...,
        Dense(config.hidden_size, 2),
    )
end

# ---------------------------------------------------------------------------
# Normalization helpers
# ---------------------------------------------------------------------------

@inline _norm(v, lo, hi) = 2 * (v - lo) / (hi - lo) - 1
@inline _denorm(v, lo, hi) = (v + 1) / 2 * (hi - lo) + lo

"""
    normalize_pinn_input(τ, T₀, A₀, η, τπ, config) -> Vector{Float64}

Map physical inputs to the [-1, 1] hypercube expected by the network.
"""
function normalize_pinn_input(
        τ::Real, T₀::Real, A₀::Real, η::Real, τπ::Real,
        config::PINNConfig,
    )
    return Float64[
        _norm(τ, config.τ_range...),
        _norm(T₀, config.T_range...),
        _norm(A₀, config.A_range...),
        _norm(η, config.η_range...),
        _norm(τπ, config.τπ_range...),
    ]
end

"""
    denormalize_pinn_output(y_norm, config) -> Vector

Convert the network output from [-1,1] back to physical [T (fm⁻¹), A].
"""
function denormalize_pinn_output(y_norm, config::PINNConfig)
    return [
        _denorm(y_norm[1], config.T_range...),
        _denorm(y_norm[2], config.A_range...),
    ]
end

# ---------------------------------------------------------------------------
# Prediction API
# ---------------------------------------------------------------------------

"""
    pinn_predict(nn, ps, st, τ, T₀, A₀, η, τπ, config) -> Vector{Float64}

Predict physical state [T(τ), A(τ)] using the trained PINN.

All inputs are in physical units; normalization is applied internally.
"""
function pinn_predict(
        nn, ps, st,
        τ::Real, T₀::Real, A₀::Real, η::Real, τπ::Real,
        config::PINNConfig,
    )
    x = normalize_pinn_input(τ, T₀, A₀, η, τπ, config)
    y_norm, _ = nn(x, ps, st)
    return denormalize_pinn_output(y_norm, config)
end
