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
    if parallel === :threads
        solutions = Vector{Any}(undef, length(initial_conditions))
        Threads.@threads for index in eachindex(initial_conditions)
            solutions[index] = solve_hydro(
                model,
                initial_conditions[index],
                tspan;
                solver = solver,
                abstol = abstol,
                reltol = reltol,
                saveat = saveat,
            )
        end
        return solutions
    else
        return [
            solve_hydro(
                model,
                u0,
                tspan;
                solver = solver,
                abstol = abstol,
                reltol = reltol,
                saveat = saveat,
            ) for u0 in initial_conditions
        ]
    end
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
