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

Wersja dla HJSW generująca warunki początkowe [T0, A0, B0] w postaci 3D SVektorów.
"""
function generate_initial_conditions(
        model::HJSWModel,
        n::Integer = 5000;
        T_range::Tuple{<:Real, <:Real} = (200.0, 1400.0),
        A_range::Tuple{<:Real, <:Real} = (-1.0, 10.0),
        B_range::Tuple{<:Real, <:Real} = (-2.0, 2.0),
        temperature_unit::Symbol = :MeV,
        seed::Integer = 5,
    )
    rng_local = Xoshiro(seed)
    ics = Vector{SVector{3, Float64}}(undef, n)

    for i in eachindex(ics)
        T0 = rand(rng_local) * (T_range[2] - T_range[1]) + T_range[1]
        T0_fm = temperature_to_fm(T0, temperature_unit)
        A0 = rand(rng_local) * (A_range[2] - A_range[1]) + A_range[1]
        B0 = rand(rng_local) * (B_range[2] - B_range[1]) + B_range[1]

        ics[i] = SVector{3, Float64}(T0_fm, A0, B0)
    end
    return ics
end
