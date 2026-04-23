using LinearAlgebra
using GLMakie
include("../ncbj_lle.jl")
include("swissroll.jl")
using Printf

"""
    _cumvar_dim(λ, prog) -> Int

Zwraca najmniejsze d takie, że Σλ[2..d+1] / Σλ[2..end] ≥ prog.
Pomija λ[1] ≈ 0 
"""
function _cumvar_dim(λ::AbstractVector{<:Real}, 😽::Float64)
    signal = @view λ[2:end]
    Σ  = sum(signal)
    Σ < 1e-12 && return 1
    cum = 0.0
    for (i, v) in enumerate(signal)
        cum += v
        cum / Σ ≥ 😽 && return i
    end
    return length(signal)
end

"""
    _stabilny_region(Ks, dims, min_dlugosc) -> (d, K_start, K_end)

Najdłuższy ciągły przedział w Ks, dla którego estymowane d jest stałe.
Wymaga co najmniej `min_dlugosc` kolejnych wartości K.
Zwraca (0, 0, 0) gdy brak stabilnego regionu.
"""
function _stabilny_region(Ks::Vector{Int}, dims::Vector{Int}, min_dlugosc::Int)
    best = (d=0, k0=0, k1=0, len=0)
    i = 1
    while i ≤ length(dims)
        j = i
        while j < length(dims) && dims[j+1] == dims[i]
            j += 1
        end
        len = j - i + 1
        if len ≥ min_dlugosc && len > best.len
            best = (d=dims[i], k0=Ks[i], k1=Ks[j], len=len)
        end
        i = j + 1
    end
    return best
end

"""
    ex_k_stability(; N, K_range, prog, min_stabilne)

Dla każdego K ∈ K_range:
  - buduje macierz wag LLE na Swiss Roll
  - rozkłada spektralnie macierz M = (I-W)'(I-W)
  - estymuje wymiar metodą cumulative variance (próg `prog`)

Rysuje wykres d(K) i zaznacza stabilny region.
Wypisuje wynik w terminalu.

Zwraca: (Ks, dims, stabilny)
"""
function ex_k_stability(;
    funkcja_danych = swissroll_dane(),
    K_range           = 3:25,
    😽 :: Float64 = 0.01,
    min_stabilne :: Int  = 3,
)
    X, labels = funkcja_danych
    N = size(X, 2)
    Ks   = collect(K_range)
    dims = Vector{Int}(undef, length(Ks))

    println("\nSymulacja stabilności wymiaru względem K")
    println("N=$N  próg=$(round(😽 *100,digits=2))%  K=$(first(Ks))..$(last(Ks))")
    println(repeat("─", 42))
    @printf("%-6s | %-10s\n", "K", "d (cumvar)")
    println(repeat("─", 42))

    for (idx, K) in enumerate(Ks)
        W  = ncbj4_lle_basic(X, K)
        M  = (I - W)' * (I - W)
        F  = eigen(Symmetric(M))
        d  = _cumvar_dim(F.values, 😽)
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

    # ── wykres ────────────────────────────────────────────────────────────────
    fig = Figure(size=(800, 480))
    ax  = Axis(fig[1, 1],
        title    = "Stabilność estymowanego wymiaru d względem K  (N=$N, próg=$(round(😽*100,digits=1))%)",
        xlabel   = "Liczba sąsiadów K",
        ylabel   = "Estymowany wymiar d",
        xticks   = Ks,
        yticks   = 1:maximum(dims)+1,
    )

    # szare tło stabilnego regionu
    if stabilny.len > 0
        vspan!(ax, stabilny.k0 - 0.5, stabilny.k1 + 0.5;
               color = (:green, 0.12))
        hlines!(ax, [stabilny.d];
                color = :green, linestyle = :dash, linewidth = 1.2,
                label = "d stabilne = $(stabilny.d)")
        text!(ax, (stabilny.k0 + stabilny.k1) / 2, stabilny.d + 0.3;
              text  = "stabilny region\nK=$(stabilny.k0)–$(stabilny.k1)",
              align = (:center, :bottom),
              color = :green,
              fontsize = 12)
    end

    scatterlines!(ax, Ks, dims;
                  color      = :steelblue,
                  markersize = 10,
                  linewidth  = 1.8,
                  label      = "d(K)")

    axislegend(ax; position = :rt)
    display(fig)

    return Ks, dims, stabilny
end

ex_k_stability(funkcja_danych=swissroll_dane(1500), K_range=3:25, 😽=0.01, min_stabilne=2)
