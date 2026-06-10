using Random
using StaticArrays

"""
Generate random initial conditions [T, A] as stack-friendly SVectors.
Temperature is always stored in internal `fm` units.
By default, points are generated using `Xoshiro(seed)` with `seed=5`.
```julia
function generate_initial_conditions(
    n::Integer;
    T_range::Tuple{<:Real,<:Real}=(400.0, 2500.0),
    A_range::Tuple{<:Real,<:Real}=(-8.0, 20.0),
    temperature_unit::Symbol=:MeV,
    seed::Integer=5,
    rng::Union{AbstractRNG,Nothing}=nothing,
)
```
"""
function generate_initial_conditions(
        n::Integer;
        T_range::Tuple{<:Real, <:Real} = (400.0, 2500.0),
        A_range::Tuple{<:Real, <:Real} = (-8.0, 20.0),
        temperature_unit::Symbol = :MeV,
        seed::Integer = 5,
        rng::Union{AbstractRNG, Nothing} = nothing,
    )
    @assert n > 0 "n must be positive."
    T_min, T_max = T_range
    A_min, A_max = A_range
    @assert T_min < T_max "Invalid T_range."
    @assert A_min < A_max "Invalid A_range."

    rng_local = isnothing(rng) ? Xoshiro(seed) : rng

    ics = Vector{SVector{2, Float64}}(undef, n)
    @inbounds for i in eachindex(ics)
        T0 = rand(rng_local) * (T_max - T_min) + T_min
        T0_fm = temperature_to_fm(T0, temperature_unit)
        A0 = rand(rng_local) * (A_max - A_min) + A_min
        ics[i] = SVector{2, Float64}(T0_fm, A0)
    end
    return ics
end

"""
    generate_initial_conditions(model::HJSWModel, n::Integer; kwargs...)

Wersja dla HJSW generująca warunki początkowe [T0, A_MIS0, A_QNM0, dA_QNM0].
"""
function generate_initial_conditions(
        model::HJSWModel,
        n::Integer = 5000;
        T_range::Tuple{<:Real, <:Real} = (200.0, 1400.0),
        A_MIS_range::Tuple{<:Real, <:Real} = (-1.0, 10.0),
        A_QNM_range::Tuple{<:Real, <:Real} = (-2.0, 2.0),
        temperature_unit::Symbol = :MeV,
        seed::Integer = 5,
        rng::Union{AbstractRNG, Nothing} = nothing,
    )
    @assert n > 0 "n must be positive."
    rng_local = isnothing(rng) ? Xoshiro(seed) : rng
    ics = Vector{SVector{4, Float64}}(undef, n)

    for i in eachindex(ics)
        T0 = rand(rng_local) * (T_range[2] - T_range[1]) + T_range[1]
        T0_fm = temperature_to_fm(T0, temperature_unit)
        A_MIS0 = rand(rng_local) * (A_MIS_range[2] - A_MIS_range[1]) + A_MIS_range[1]
        A_QNM0 = rand(rng_local) * (A_QNM_range[2] - A_QNM_range[1]) + A_QNM_range[1]
        dA_QNM0 = rand(rng_local) * 0.1

        ics[i] = SVector{4, Float64}(T0_fm, A_MIS0, A_QNM0, dA_QNM0)
    end
    return ics
end
