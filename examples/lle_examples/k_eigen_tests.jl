using LinearAlgebra
using GLMakie
using Printf

include("../ncbj_lle.jl")
include("swissroll.jl")
include("2d_examples_data.jl")

function _eigengap_dim(λ::AbstractVector{<:Real}; max_d::Int = 10)
    limit = min(length(λ), max_d + 2)
    λ_robocze = @view λ[2:limit]
    przyrosty = diff(λ_robocze)
    return argmax(przyrosty)
end

function _stabilny_region(Ks::Vector{Int}, dims::Vector{Int}, min_dlugosc::Int)
    best = (d = 0, k0 = 0, k1 = 0, len = 0)
    i = 1
    while i ≤ length(dims)
        j = i
        while j < length(dims) && dims[j + 1] == dims[i]
            j += 1
        end
        len = j - i + 1
        if len ≥ min_dlugosc && len > best.len
            best = (d = dims[i], k0 = Ks[i], k1 = Ks[j], len = len)
        end
        i = j + 1
    end
    return best
end

function ex_k_stability(;
        data,
        K_range = 3:25,
        max_d::Int = 10,
        min_stabilne::Int = 3,
    )
    X = data
    N = size(X, 2)
    Ks = collect(K_range)
    dims = Vector{Int}(undef, length(Ks))

    println("\nSymulacja stabilności wymiaru względem K")
    println("N=$N  K=$(first(Ks))..$(last(Ks))  (Metoda Eigengap)")
    println(repeat("─", 42))
    @printf("%-6s | %-10s\n", "K", "d (eigengap)")
    println(repeat("─", 42))

    for (idx, K) in enumerate(Ks)
        W = ncbj4_lle_basic(X, K)
        M = (I - W)' * (I - W)
        F = eigen(Symmetric(M))

        # Użycie nowej metody wyznaczania wymiaru
        d = _eigengap_dim(F.values, max_d = max_d)

        dims[idx] = d
        @printf("%-6d | %-10d\n", K, d)
    end

    stabilny = _stabilny_region(Ks, dims, min_stabilne)

    println(repeat("─", 42))
    if stabilny.len > 0
        println("Stabilny region: K=$(stabilny.k0)..$(stabilny.k1)  →  d = $(stabilny.d)  ($(stabilny.len) kolejnych K)")
    else
        println("Brak stabilnego regionu przy min_stabilne=$min_stabilne")
    end

    fig = Figure(size = (800, 480))
    ax = Axis(
        fig[1, 1],
        title = "Stabilność estymowanego wymiaru d (Eigengap) względem K (N=$N)",
        xlabel = "Liczba sąsiadów K",
        ylabel = "Estymowany wymiar d",
        xticks = Ks,
        yticks = 1:(maximum(dims) + 1),
    )

    scatterlines!(
        ax, Ks, dims;
        color = :steelblue,
        markersize = 10,
        linewidth = 1.8,
        label = "d(K)"
    )

    axislegend(ax; position = :rt)
    display(fig)

    return Ks, dims, stabilny
end
