using Lux
using ForwardDiff
using LinearAlgebra
using Statistics

"""
    normalize_generic(X::AbstractMatrix{<:Real})

Pomocnik do dynamicznego mapowania dowolnej liczby kolumn danych na zakres [-1, 1].
Zwraca znormalizowane dane oraz wektory minimum i maksimum dla każdego wymiaru.
"""
function normalize_generic(X::AbstractMatrix{<:Real})
    mn = minimum(X, dims = 1)
    mx = maximum(X, dims = 1)
    r = mx .- mn
    r[r .== 0.0] .= 1.0
    X_norm = 2.0 .* (X .- mn) ./ r .- 1.0
    return X_norm, vec(mn), vec(mx)
end

"""
    build_generic_network(input_dim::Int, output_dim::Int; hidden_layers=3, hidden_size=64)

Buduje uniwersalną sieć MLP (Chain) dopasowaną do liczby kolumn w Twoich danych.
"""
function build_generic_network(input_dim::Int, output_dim::Int; hidden_layers = 3, hidden_size = 64)
    hidden = ntuple(_ -> Dense(hidden_size, hidden_size, tanh), hidden_layers - 1)
    return Chain(Dense(input_dim, hidden_size, tanh), hidden..., Dense(hidden_size, output_dim))
end

"""
    compute_data_jacobian(nn, ps, st, τ, x0, mn_in, mx_in, mn_out, mx_out) -> Matrix{Float64}

Oblicza Jakobian ewolucji bezpośrednio z wag sieci wytrenowanej na danych eksperymentalnych.
Automatycznie adaptuje się do dowolnego wymiaru N.

# Argumenty:
- `nn`, `ps`, `st`: sieć Lux, parametry i stan
- `τ`: proper time (fizyczny czas)
- `x0`: stan początkowy w chwili początkowej (fizyczny, bez czasu tau)
- `mn_in`, `mx_in`: minima i maksima zmiennych wejściowych [tau, x0...]
- `mn_out`, `mx_out`: minima i maksima zmiennych wyjściowych
"""
function compute_data_jacobian(
        nn, ps, st,
        τ::Real, x0::AbstractVector{<:Real},
        mn_in::AbstractVector{<:Real}, mx_in::AbstractVector{<:Real},
        mn_out::AbstractVector{<:Real}, mx_out::AbstractVector{<:Real},
    )
    # 1. Normalizacja wejścia
    τ_n = 2.0 * (τ - mn_in[1]) / (mx_in[1] - mn_in[1]) - 1.0
    x0_n = 2.0 .* (x0 .- mn_in[2:end]) ./ (mx_in[2:end] - mn_in[2:end]) .- 1.0

    # 2. Funkcja różniczkowalna przez ForwardDiff
    f = x0_vec_n -> begin
        x_input = vcat(τ_n, x0_vec_n)
        y_norm, _ = nn(x_input, ps, st)
        return y_norm
    end

    # 3. Automatyczny Jakobian w przestrzeni znormalizowanej
    J_norm = ForwardDiff.jacobian(f, x0_n)

    # 4. Skalowanie powrotne do jednostek fizycznych
    scale_out = (mx_out .- mn_out) ./ 2.0
    inv_scale_in = 2.0 ./ (mx_in[2:end] .- mn_in[2:end])

    J_phys = J_norm .* (scale_out .* inv_scale_in')
    return J_phys
end
