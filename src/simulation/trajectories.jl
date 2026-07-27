using DifferentialEquations

"""
    generate_trajectories(model, initial_conditions, tspan; saveat=nothing, parallel=:serial)

Generuje wektor rozwiązań równań hydrodynamiki dla zadanych warunków początkowych.
"""
function generate_trajectories(
        model::AbstractHydroModel,
        initial_conditions::AbstractVector{<:AbstractVector{<:Real}},
        tspan::Tuple{<:Real, <:Real};
        saveat = nothing,
        parallel::Symbol = :serial,
        abstol::Real = 1.0e-8,
        reltol::Real = 1.0e-8,
        solver = Rodas5(),
    )
    # 1. Bazowy problem różniczkowy dla 1. warunku początkowego
    Tstate = promote_type(eltype(initial_conditions[1]), typeof(tspan[1]), typeof(tspan[2]), Float64)
    N = length(initial_conditions[1])
    u0_first = SVector{N, Tstate}(initial_conditions[1]...)
    base_prob = ODEProblem(rhs, u0_first, (Tstate(tspan[1]), Tstate(tspan[2])), model)

    # 2. Definicja problemu zespołowego (kompatybilna ze wszystkimi wersjami SciMLBase)
    function prob_func(prob, arg2, args...)
        i = isempty(args) ? (arg2 isa Integer ? arg2 : arg2.sim_id) : arg2
        return remake(prob, u0 = initial_conditions[i])
    end
    ensemble_prob = EnsembleProblem(base_prob; prob_func = prob_func)

    # 3. Dobór algorytmu równoległości
    ensemble_alg = parallel === :threads ? EnsembleThreads() :
                   parallel === :split_threads ? EnsembleSplitThreads() : EnsembleSerial()

    # 4. Rozwiązanie zespołu trajektorii
    extra_kwargs = isnothing(saveat) ? (;) : (; saveat)
    ens_sol = solve(ensemble_prob, solver, ensemble_alg;
                    trajectories = length(initial_conditions),
                    abstol = abstol, reltol = reltol, extra_kwargs...)

    return ens_sol.u
end

function build_dataset(solutions::AbstractVector; temperature_unit::Symbol = :fm)
    state_dim = length(solutions[1].u[1])
    n_cols = 1 + state_dim
    n_rows = sum(length(sol.t) for sol in solutions)
    data = Matrix{Float64}(undef, n_rows, n_cols)

    row_idx = 1
    @inbounds for sol in solutions
        for i in eachindex(sol.t)
            tau = sol.t[i]
            all_finite = isfinite(tau) && all(isfinite, sol.u[i])

            if all_finite
                data[row_idx, 1] = tau

                T = sol.u[i][1]
                if temperature_unit === :fm
                    # standard fm units
                elseif temperature_unit === :MeV
                    T = T * MEV_PER_FM
                else
                    throw(ArgumentError("Unsupported temperature_unit: $temperature_unit. Expected :fm or :MeV."))
                end
                data[row_idx, 2] = T

                for col in 3:n_cols
                    data[row_idx, col] = sol.u[i][col - 1]
                end

                row_idx += 1
            end
        end
    end

    return data[1:(row_idx - 1), :]
end
