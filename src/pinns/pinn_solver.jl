"""
This file provides:
- `predict_trajectories_pinn`  – predict solution families, compatible with `build_dataset`
- `build_pinn_dataset`          – convenience wrapper
- `compare_pinn_ode`            – point-wise accuracy diagnostics
- `pinn_attractor_analysis`     – run PCA on PINN predictions to track dimensionality
"""

"""
    predict_trajectories_pinn(result, initial_conditions, model, tspan; saveat)

Predict a family of trajectories using the trained PINN.

Returns a `Vector` of `NamedTuple` objects `(t, u)` that are fully
compatible with the existing `build_dataset` function.

# Arguments
- `result`:             trained `PINNResult`
- `initial_conditions`: `Vector` of SVectors `[T₀, A₀]`  (fm units)
- `model`:              `AbstractHydroModel` (used for η/s and τ_π)
- `tspan`:              `(τ_start, τ_end)` in fm/c
- `saveat`:             step size or vector of save times (default 0.01)
"""
function predict_trajectories_pinn(
        result::PINNResult,
        initial_conditions::AbstractVector,
        model::AbstractHydroModel,
        tspan::Tuple{<:Real, <:Real};
        saveat::Union{Real, AbstractVector{<:Real}, Nothing} = 0.01,
    )
    τ_grid = if saveat isa Real
        collect(range(tspan[1], tspan[2]; step = saveat))
    elseif saveat isa AbstractVector
        collect(Float64, saveat)
    else
        range(tspan[1], tspan[2]; length = 100) |> collect
    end

    η = model.params.eta_over_s
    τπ = model.params.tau_pi

    nn, ps, st, cfg = result.network, result.ps, result.st, result.config

    return map(initial_conditions) do ic
        T₀, A₀ = Float64(ic[1]), Float64(ic[2])
        u_vals = [pinn_predict(nn, ps, st, τ, T₀, A₀, η, τπ, cfg) for τ in τ_grid]
        (t = τ_grid, u = u_vals)
    end
end


"""
    build_pinn_dataset(result, initial_conditions, model, tspan; saveat, temperature_unit)

Build a `[n_rows × 3]` dataset `[τ, T, A]` from PINN predictions.

Drop-in replacement for `build_dataset` that uses the trained PINN instead
of the ODE solver.  The result can be passed directly to `run_pca_per_time`,
`run_evolution_pca_workflow`, `estimate_dimension`, etc.
"""
function build_pinn_dataset(
        result::PINNResult,
        initial_conditions::AbstractVector,
        model::AbstractHydroModel,
        tspan::Tuple{<:Real, <:Real};
        saveat::Union{Real, AbstractVector{<:Real}, Nothing} = 0.01,
        temperature_unit::Symbol = :fm,
    )
    pinn_sols = predict_trajectories_pinn(
        result, initial_conditions, model, tspan; saveat = saveat,
    )
    return build_dataset(pinn_sols; temperature_unit = temperature_unit)
end


"""
    compare_pinn_ode(result, model, ic, tspan; saveat) -> NamedTuple

Compare PINN predictions to the reference ODE solver for one trajectory.

Returns a `NamedTuple` with:
- `tau`     : time grid
- `T_ode`, `A_ode`   : reference solution
- `T_pinn`, `A_pinn` : PINN prediction
- `T_error`, `A_error` : absolute point-wise errors

# Example
```julia
result = train_pinn(model, config; n_epochs = 3000)
ic = [2.0, 0.5]   # [T₀ (fm⁻¹), A₀]
cmp = compare_pinn_ode(result, model, ic, (0.22, 1.2))
```
"""
function compare_pinn_ode(
        result::PINNResult,
        model::AbstractHydroModel,
        ic::AbstractVector,
        tspan::Tuple{<:Real, <:Real};
        saveat::Union{Real, AbstractVector{<:Real}, Nothing} = 0.01,
    )
    # Reference ODE
    ode_sol = solve_hydro(model, ic, tspan; saveat = saveat)
    τ_grid = ode_sol.t

    T₀, A₀ = Float64(ic[1]), Float64(ic[2])
    η = model.params.eta_over_s
    τπ = model.params.tau_pi
    nn, ps, st, cfg = result.network, result.ps, result.st, result.config

    # PINN predictions at the same time points
    pinn_vals = [pinn_predict(nn, ps, st, τ, T₀, A₀, η, τπ, cfg) for τ in τ_grid]

    T_ode = [u[1] for u in ode_sol.u]
    A_ode = [u[2] for u in ode_sol.u]
    T_pinn = [v[1] for v in pinn_vals]
    A_pinn = [v[2] for v in pinn_vals]

    return (
        tau = τ_grid,
        T_ode = T_ode, A_ode = A_ode,
        T_pinn = T_pinn, A_pinn = A_pinn,
        T_error = abs.(T_pinn .- T_ode),
        A_error = abs.(A_pinn .- A_ode),
    )
end


"""
    pinn_attractor_analysis(result, model, tspan; n_points, saveat, n_components, kwargs...)
    -> NamedTuple

Full attractor-analysis workflow driven by PINN predictions.

Generates `n_points` random initial conditions (uniformly over the config
ranges), predicts trajectories with the PINN, then runs `run_pca_per_time`
to track the evolution of explained variance ratio (EVR) over proper time.

When EVR of the first PC approaches 1, the system has collapsed onto a
low-dimensional attractor.

# Returns
Named tuple with:
- `initial_conditions`
- `pinn_solutions`
- `dataset`
- `pca_over_time`  (same layout as `run_evolution_pca_workflow`)

# Example
```julia
result = train_pinn(BRSSSModel(), PINNConfig(); n_epochs = 5000)
ana    = pinn_attractor_analysis(result, BRSSSModel(), (0.22, 1.2))
# ana.pca_over_time.explained_variance_ratio[:,1]  → PC1 share vs τ
```
"""
function pinn_attractor_analysis(
        result::PINNResult,
        model::AbstractHydroModel,
        tspan::Tuple{<:Real, <:Real};
        n_points::Int = 500,
        saveat = 0.01,
        n_components::Int = 2,
        seed::Int = 7,
        feature_cols::AbstractVector{<:Integer} = [2, 3],
        temperature_unit::Symbol = :fm,
    )
    cfg = result.config
    rng = MersenneTwister(seed)

    _rnd(lo, hi) = lo + rand(rng) * (hi - lo)

    # Sample ICs from config ranges (T in fm⁻¹, A dimensionless)
    ics = [
        Float64[_rnd(cfg.T_range...), _rnd(cfg.A_range...)]
            for _ in 1:n_points
    ]

    pinn_sols = predict_trajectories_pinn(
        result, ics, model, tspan; saveat = saveat,
    )

    dataset = build_dataset(pinn_sols; temperature_unit = temperature_unit)

    pca_over_time = run_pca_per_time(
        dataset;
        n_components = n_components,
        feature_cols = feature_cols,
    )

    return (
        initial_conditions = ics,
        pinn_solutions = pinn_sols,
        dataset = dataset,
        pca_over_time = pca_over_time,
    )
end
