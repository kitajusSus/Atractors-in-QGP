using Printf
using LinearAlgebra
using DelimitedFiles
include("../../examples/ncbj_lle.jl")

function compare_eigen_jlmatlab(; K = 12, indeksy = 1:10)
    X_matlab = readdlm(joinpath(@__DIR__, "X_matlab.csv"))
    E_matlab = vec(readdlm(joinpath(@__DIR__, "E_matlab.csv")))

    W_julia = ncbj4_lle_basic(X_matlab, K)
    M_julia = (I - W_julia)' * (I - W_julia)
    F_julia = eigen(Symmetric(M_julia))

    val_matlab = sort(E_matlab)
    # roznie bywa z tym matlabem wiec bierze tyle ile matlab/octave zwraca
    limit = min(length(val_matlab), length(indeksy))
    idx = indeksy[1:limit]

    val_m = val_matlab[idx]
    val_j = F_julia.values[idx]

    norm_m = val_m ./ val_m[2]
    norm_j = val_j ./ val_j[2]

    println("\nPorównanie dla K = $K")
    @printf("%-8s | %-16s | %-16s | %-16s | %-12s | %-12s\n", "Indeks", "Matlab (surowe)", "Julia (surowe)", "Różnica (sur.)", "Matlab (nor)", "Julia (nor)")
    println(repeat("=", 70))

    for i in 1:length(idx)
        diff_surowe = abs(val_m[i] - val_j[i])
        @printf(
            "%-8d | %-16.6e | %-16.6e | %-16.6e | %-12.4f | %-12.4f\n",
            idx[i], val_m[i], val_j[i], diff_surowe, norm_m[i], norm_j[i]
        )
    end

    return norm_m, norm_j
end

compare_eigenvalues(K = 12, indeksy = 1:4)
