using Printf
using LinearAlgebra
using GLMakie

include("../ncbj_lle.jl")
include("swissroll.jl")

function lambda_analiza(; X, K=12, indeksy=1:10)
    
    
    W = ncbj4_lle_basic(X, K)
    M = (I - W)' * (I - W)
    F = eigen(Symmetric(M))
    
    limit = min(length(F.values), maximum(indeksy))
    idx = indeksy[indeksy .<= limit]
    
    val_j = F.values[idx]
    norm_j = val_j
    norm_j[2:end] = val_j[2:end] ./ val_j[2]
    
    # println("\nAnaliza widma dla K = $K")
    # @printf("%-8s | %-18s | %-18s\n", "Indeks", "Wartość surowa", "Wartość znormalizowana")
    # println(repeat("=", 52))
    #
    # for i in 1:length(idx)
    #     @printf("%-8d | %-18.6e | %-18.4f\n", idx[i], val_j[i], norm_j[i])
    # end
    
    return norm_j
end





#
# t candidates are:
#   lambda_analiza(; X, K, indeksy)
#    @ Main ~/github/Atractors-in-QGP/examples/lle_examples/k_variance_lambda.jl:8
#
# Stacktrace:
#  [1] top-level scope
#    @ REPL[13]:1
#
# julia> lambda_analiza(X= ξ, K=25, indeksy=1:10)
#
# Analiza widma dla K = 25
# Indeks   | Wartość surowa     | Wartość znormalizowana
# ====================================================
# 1        | -3.887918e-17      | -0.0000
# 2        | 1.104921e-09       | 1.0000
# 3        | 1.184392e-07       | 107.1925
# 4        | 1.413256e-07       | 127.9057
# 5        | 5.656513e-07       | 511.9382
# 6        | 1.166524e-06       | 1055.7534
# 7        | 3.281763e-06       | 2970.1339
# 8        | 6.405369e-06       | 5797.1286
# 9        | 1.057633e-05       | 9572.0279
# 10       | 1.765323e-05       | 15976.9177
# ([-3.518729001240954e-8, 1.0, 107.19247249091299, 127.90565364000673, 511.9381998913114, 1055.7534229182668, 2970.1338659767052, 5797.128575266762, 9572.027877470493, 15976.91765595359], [-3.887917553222838e-17, 1.104920996147098e-9, 1.1843921348412998e-7, 1.4132564223276193e-7, 5.6565126578966e-7, 1.1665241237365599e-6, 3.2817632698852127e-6, 6.4053690801765574e-6, 1.057633457752249e-5, 1.7653231771776398e-5])
#
# julia> 𝐗 =lambda_analiza(X= ξ, K=25, indeksy=1:10)
# 10-element Vector{Float64}:
#     -3.518729001240954e-8
#      1.0
#    107.19247249091299
#    127.90565364000673
#    511.9381998913114
#   1055.7534229182668
#   2970.1338659767052
#   5797.128575266762
#   9572.027877470493
#  15976.91765595359
#
# julia> \ma
#
# ♂ \male      ✠ \maltese
# ↧ \mapsdown  ↤ \mapsfrom
# ↦ \mapsto    ↥ \mapsup
# ♂ \mars
# julia> ♌ = zeros(Float64, 40,10)
# 40×10 Matrix{Float64}:
#  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
#  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
#  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
#  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
#  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
#  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
#  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
#  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
#  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
#  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
#  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
#  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
#  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
#  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
#  ⋮                        ⋮
#  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
#  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
#  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
#  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
#  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
#  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
#  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
#  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
#  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
#  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
#  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
#  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
#  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
#
# julia> for i in 1:45
#            ♌[i,:] = lambda_analiza()
# lambda_analiza(; X, K, indeksy) @ Main ~/github/Atractors-in-QGP/examples/lle_examples/k_variance_lambda.jl:8
# julia> for i in 1:45
#            ♌[i,:] = lambda_analiza()
# lambda_analiza(; X, K, indeksy) @ Main ~/github/Atractors-in-QGP/examples/lle_examples/k_variance_lambda.jl:8
# julia> for i in 1:45
#            ♌[i,:] = lambda_analiza(X = ξ, K=i, indeksy=1:10)
#        end
# """
#
