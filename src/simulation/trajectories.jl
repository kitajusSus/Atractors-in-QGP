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
        parallel = nothing,
        abstol::Real = 1.0e-8,
        reltol::Real = 1.0e-8,
        solver = Rodas5(),
    )
    first_initial_state = SVector{length(initial_conditions[1]), Float64}(initial_conditions[1]...)
    base_problem = ODEProblem(rhs, first_initial_state, tspan, model)

    function problem_generator_function(problem, second_argument, remaining_arguments...)
        trajectory_index = isempty(remaining_arguments) ? (second_argument isa Integer ? second_argument : second_argument.sim_id) : second_argument
        return remake(problem, u0 = initial_conditions[trajectory_index])
    end

    ensemble_problem = EnsembleProblem(base_problem; prob_func = problem_generator_function)

    extra_keyword_arguments = isnothing(saveat) ? (;) : (; saveat)
    ensemble_solution = solve(
        ensemble_problem,
        solver,
        EnsembleThreads();
        trajectories = length(initial_conditions),
        abstol = abstol,
        reltol = reltol,
        extra_keyword_arguments...
    )

    return ensemble_solution.u
end

function build_dataset(solutions::AbstractVector; temperature_unit::Symbol = :fm, has_temperature::Bool = true)
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

                if has_temperature
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
                else
                    for col in 2:n_cols
                        data[row_idx, col] = sol.u[i][col - 1]
                    end
                end

                row_idx += 1
            end
        end
    end

    return data[1:(row_idx - 1), :]
end

