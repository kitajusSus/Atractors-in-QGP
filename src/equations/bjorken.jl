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

"""
RHS dla przepływu Bjorkena w teorii HJSW.
Układ współrzędnych: u = [T, A_MIS, A_QNM, dA_QNM]
"""
function rhs(u::AbstractVector{<:Real}, model::HJSWModel, τ::Real)
    # 𝚷_MIS = A_MIS * T^4
    # 𝚷_QNM = A_QNM * T^4
    #
    T, A_MIS, A_QNM, dA_QNM = u[1], u[2], u[3], u[4]
    p = model.params
    # 𝚷_total = 𝚷_MIS + 𝚷_QNM

    A_total = A_MIS + A_QNM

    dT = (T / τ) * (-1.0 / 3.0 + A_total / 6.0)
    # navier stokes
    A_NS = 8 / (3 * τ * T) * p.eta_over_s
    dA_MIS = (A_NS - A_MIS) / (p.tau_pi / T) - (4 * A_MIS) / (3 * τ)
    # mod holograficzny
    ω2 = model.omega_R^2 + model.omega_I^2

    d2A_QNM = (
        - (2 * model.omega_I * T + 4 / (3 * τ) - dT / T) * dA_QNM
            - (ω2 * T^2 + (4 * model.omega_I * T) / (3 * τ) - 2 / (9 * τ^2) - (2 * dT) / (3 * τ * T)) * A_QNM
    )

    return SVector(dT, dA_MIS, dA_QNM, d2A_QNM)
end
