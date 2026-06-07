using DifferentialEquations
using ProgressMeter

function generate_trajectories(
        model::AbstractHydroModel,
        initial_conditions::AbstractVector{<:AbstractVector{<:Real}},
        tspan::Tuple{<:Real, <:Real};
        saveat = nothing,
        parallel::Symbol = :serial,
    )
    @assert !isempty(initial_conditions) "At least one initial condition is required."

    first_solution = solve_hydro(model, initial_conditions[1], tspan; saveat = saveat)
    solType = typeof(first_solution)
    solutions = Vector{solType}(undef, length(initial_conditions))
    solutions[1] = first_solution

    if length(initial_conditions) == 1
        return solutions
    end

    base_problem = first_solution.prob

    function prob_func(prob, arg2, args...)
        i = if isempty(args)
            ctx = arg2
            if hasproperty(ctx, :sim_id)
                ctx.sim_id
            else
                error("Nieznana struktura EnsembleContext. Dostępne pola: $(propertynames(ctx))")
            end
        else
            arg2
        end
        return remake(prob, u0 = initial_conditions[i])
    end

    pasek = Progress(length(initial_conditions), 1, "Obliczanie trajektorii: ")

    function output_func(sol, arg2, args...)
        next!(pasek)
        return (sol, false)
    end

    ensemble_problem = EnsembleProblem(base_problem; prob_func = prob_func, output_func = output_func)
    ensemble_alg = parallel === :threads ? EnsembleThreads() : EnsembleSerial()

    solve_kwargs = (
        trajectories = length(initial_conditions),
        abstol = 1.0e-8,
        reltol = 1.0e-8,
    )

    ensemble_solution = if isnothing(saveat)
        solve(ensemble_problem, Rodas5(), ensemble_alg; solve_kwargs...)
    else
        solve(ensemble_problem, Rodas5(), ensemble_alg; solve_kwargs..., saveat = saveat)
    end

    @inbounds for i in eachindex(solutions)
        solutions[i] = ensemble_solution.u[i]
    end

    return solutions
end

function build_dataset(solutions::AbstractVector; temperature_unit::Symbol = :fm)
    n_rows = sum(length(sol.t) for sol in solutions)
    data = Matrix{Float64}(undef, n_rows, 3)
    hbarc = 197.3269804

    row_idx = 1
    @inbounds for sol in solutions
        for i in eachindex(sol.t)
            tau = sol.t[i]
            T = sol.u[i][1]
            A = sol.u[i][2]

            if temperature_unit === :MeV
                T = T * hbarc
            end

            data[row_idx, 1] = tau
            data[row_idx, 2] = T
            data[row_idx, 3] = A

            row_idx += 1
        end
    end

    return data
end
