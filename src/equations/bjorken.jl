using StaticArrays

"""
In-place RHS for Bjorken flow in BRSSS-like variables u = [T, A].
defined as:


```julia
function rhs(u::AbstractVector{<:Real}, model::AbstractHydroModel, τ::Real)
end
```

Compute the right-hand side of the Bjorken flow equations for
a given state `u`, model parameters, and proper time `τ`.
```julia
function rhs!(du, u, model, τ)
    # Compute dT and dA based on the equations of motion
    du[1] = dT
    du[2] = dA
end
```

returns:

```julia
    return SVector(dT, dA)
```
"""
function rhs(u::AbstractVector{<:Real}, model::AbstractHydroModel, τ::Real)
    T, A = u[1], u[2]

    p = model.params

    dT = (T / τ) * (-one(T) / 3 + A / 18)
    term_T = τ * T * (A + (p.lambda1 / (12 * p.eta_over_s)) * A^2)

    term_A2 = (2 / 9) * p.tau_pi * A^2
    dA = (1 / (p.tau_pi * τ)) * (8 * p.eta_over_s - term_T - term_A2)
    return SVector(dT, dA)
end

function rhs(u::AbstractVector{<:Real}, model::HJSWModel, t::Real)
    T, A, B = u[1], u[2], u[3]

    Ω_I = model.omega_I
    Ω_R = model.omega_R
    C_η = model.params.eta_over_s

    Ω2 = Ω_I^2 + Ω_R^2

    dT = (1.0 / 18.0) * (A - 6.0) * T / t

    dA = B / t

    term1 = A * (- (11.0 * B) / 18.0 - t^2 * T^2 * Ω2)
    term2 = - (4.0 / 27.0) * (A^2) * (3.0 * t * Ω_I * T - 1.0)
    term3 = - (1.0 / 27.0) * (A^3)
    term4 = B * ((2.0 / 3.0) - 2.0 * t * Ω_I * T)
    term5 = 8.0 * C_η * t * T * Ω2

    dB = (term1 + term2 + term3 + term4 + term5) / t

    return SVector(dT, dA, dB)
end
