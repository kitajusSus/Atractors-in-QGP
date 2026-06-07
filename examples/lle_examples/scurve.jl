using GLMakie
using LinearAlgebra
const dx = 1.0e-3
include("../ncbj_lle.jl")
# skopiowane z examples/ncbj_lle.jl
include("plots_lle.jl")

function scurve_dane(N)
    θ = LinRange(-1.5 * π, 1.5 * π, N)
    x = sin.(θ) .+ 0.1 .* randn(Float64, N)
    y = LinRange(0.0, 5.0, N) .+ 0.01 .* randn(Float64, N)
    z = sign.(θ) .+ 0.1 .* randn(Float64, N)
    X = [x y z]'
    return X, θ
end
function ex_scurve(; N = 2000, K = 12)
    X, labels = scurve_dane(N)
    wagi = ncbj4_lle_basic(X, K)
    Y = ncbj5_nowy_manifold(wagi)
    fig = plot_examples_lle(X, Y, labels)
    # display(fig)
    return fig
end
