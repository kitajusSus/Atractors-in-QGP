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

function build_dataset(
        solutions::AbstractVector;
        temperature_unit::Symbol = :fm,
        has_temperature::Bool = true,
        include_tau::Bool = false,
        include_w::Bool = false,
        include_temperature::Bool = true,
    )
    state_dim = length(solutions[1].u[1])
    has_extra = has_temperature && (include_tau || include_w)
    
    n_cols = if has_extra
        include_temperature ? (2 + state_dim) : (1 + state_dim)
    elseif has_temperature
        include_temperature ? (1 + state_dim) : state_dim
    else
        1 + state_dim
    end

    n_rows = sum(length(sol.t) for sol in solutions)
    data = Matrix{Float64}(undef, n_rows, n_cols)

    row_idx = 1
    @inbounds for sol in solutions
        for i in eachindex(sol.t)
            t_var = sol.t[i]
            all_finite = isfinite(t_var) && all(isfinite, sol.u[i])

            if all_finite
                if has_extra
                    # Determine tau and w
                    T_fm = sol.u[i][1]
                    if include_tau
                        tau = t_var / T_fm
                        w = t_var
                    else
                        tau = t_var
                        w = tau * T_fm
                    end
                    T = temperature_unit === :MeV ? T_fm * MEV_PER_FM : T_fm

                    if include_temperature
                        data[row_idx, 1] = tau
                        data[row_idx, 2] = T
                        data[row_idx, 3] = w
                        for col in 4:n_cols
                            data[row_idx, col] = sol.u[i][col - 2]
                        end
                    else
                        data[row_idx, 1] = tau
                        data[row_idx, 2] = w
                        for col in 3:n_cols
                            data[row_idx, col] = sol.u[i][col - 1]
                        end
                    end
                elseif has_temperature
                    data[row_idx, 1] = t_var
                    T = sol.u[i][1]
                    if temperature_unit === :fm
                        # standard fm units
                    elseif temperature_unit === :MeV
                        T = T * MEV_PER_FM
                    else
                        throw(ArgumentError("Unsupported temperature_unit: $temperature_unit. Expected :fm or :MeV."))
                    end
                    
                    if include_temperature
                        data[row_idx, 2] = T
                        for col in 3:n_cols
                            data[row_idx, col] = sol.u[i][col - 1]
                        end
                    else
                        for col in 2:n_cols
                            data[row_idx, col] = sol.u[i][col]
                        end
                    end
                else
                    data[row_idx, 1] = t_var
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

