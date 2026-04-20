using LinearAlgebra
using Statistics
using CairoMakie
using SparseArrays

include("../ncbj_lle.jl")
include("plots_lle.jl")

function kula_3d_w_15d(N)
    T = randn(Float64, 3, N)
    X = zeros(Float64, 15, N)
    
    X[1:3, :] .= T
    X[4, :] .= sin.(T[1, :])
    X[5, :] .= T[2, :] .^ 2
    X[6, :] .= exp.(.-abs.(T[3, :]))
    
    A = randn(Float64, 9, 3)
    X[7:15, :] .= A * T
    
    X .+= 0.05 .* randn(Float64, 15, N)
    labels = T[1, :]
    
    return X, labels
end

function hiper_rozmaitosc_dane(N)
    T = randn(Float64, 4, N)
    X = zeros(Float64, 20, N)
    
    X[1:4, :] .= T
    X[5, :] .= sin.(T[1, :])
    X[6, :] .= cos.(T[2, :])
    X[7, :] .= T[3, :] .^ 2
    X[8, :] .= T[4, :] .^ 3
    X[9, :] .= T[1, :] .* T[2, :]
    X[10, :] .= exp.(.- abs.(T[3, :]))
    
    A = randn(Float64, 10, 4)
    X[11:20, :] .= A * T
    
    X .+= 0.05 .* randn(Float64, 20, N)
    labels = T[1, :]
    
    return X, labels
end

function buduj_macierz_M(W)
    I_W = I - W
    return I_W' * I_W
end

function oblicz_wartosci_wlasne(X, K, limit_wartosci)
    W = ncbj4_lle_basic(X, K)
    M = buduj_macierz_M(W)
    
    F = eigen(Symmetric(Matrix(M)))
    
    dostepne = min(limit_wartosci, length(F.values))
    wartosci = zeros(Float64, limit_wartosci)
    wartosci[1:dostepne] .= F.values[1:dostepne]
    
    return wartosci
end

function skanuj_wymiary_K(X, K_zakres, limit_wartosci)
    macierz_wartosci = zeros(Float64, length(K_zakres), limit_wartosci)
    
    for (i, K) in enumerate(K_zakres)
        macierz_wartosci[i, :] .= oblicz_wartosci_wlasne(X, K, limit_wartosci)
    end
    
    return macierz_wartosci
end

function oblicz_wariancje(macierz_wartosci)
    liczba_wartosci = size(macierz_wartosci, 2)
    wariancje = zeros(Float64, liczba_wartosci)
    
    for j in 1:liczba_wartosci
        kolumna = macierz_wartosci[:, j]
        wariancje[j] = var(kolumna)
    end
    
    return wariancje
end
function estymuj_wymiar_danych(wariancje::Vector{Float64})
    n = length(wariancje)
    
    if n < 3
        return 1
    end
    
    najlepszy_indeks = 2
    maks_wariancja_miedzyklasowa = -1.0
    
    for k in 2:(n-1)
        klasa_sygnalu = view(wariancje, 2:k)
        klasa_szumu = view(wariancje, (k+1):n)
        
        waga_sygnalu = length(klasa_sygnalu) / (n - 1)
        waga_szumu = length(klasa_szumu) / (n - 1)
        
        # Zapobiegamy błędom przy pustych klasach
        srednia_sygnalu = isempty(klasa_sygnalu) ? 0.0 : mean(klasa_sygnalu)
        srednia_szumu = isempty(klasa_szumu) ? 0.0 : mean(klasa_szumu)
        srednia_calkowita = mean(view(wariancje, 2:n))
        
        var_miedzy = waga_sygnalu * (srednia_sygnalu - srednia_calkowita)^2 + 
                     waga_szumu * (srednia_szumu - srednia_calkowita)^2
                     
        if var_miedzy > maks_wariancja_miedzyklasowa
            maks_wariancja_miedzyklasowa = var_miedzy
            najlepszy_indeks = k
        end
    end
    
    # Odejmujemy 1, ponieważ pierwszy indeks był zerowy i nie liczy się do wymiaru
    wymiar = najlepszy_indeks - 1
    
    return max(1, wymiar)
end
function rysuj_wariancje(wariancje, wykryty_wymiar, nazwa_zbioru)
    set_publication_theme_large()
    ile_wartosci = length(wariancje)
    fig = Figure(size=(800, 600))
    ax = Axis(fig[1, 1], 
        title="Wariancja λ dla $nazwa_zbioru (Wykryty wymiar: $wykryty_wymiar)", 
        xlabel="Indeks wartosci wlasnej", 
        ylabel="Wariancja"
    )
    scatterlines!(ax, 1:ile_wartosci, wariancje, markersize=12, color=:red)
    display(fig)
    return fig
end

function analizuj_zbior_danych(X_dane, nazwa_zbioru; K_zakres=1:1:20, limit_wartosci=20)
    println("Rozpoczynam analize zbioru: $nazwa_zbioru")
    
    macierz_wartosci = skanuj_wymiary_K(X_dane, K_zakres, limit_wartosci)
    wariancje = oblicz_wariancje(macierz_wartosci)
    wykryty_wymiar = estymuj_wymiar_danych(wariancje)
    
    println("Zakonczono. Wykryty wymiar danych to: ", wykryty_wymiar, "\n")
    
    rysuj_wariancje(wariancje, wykryty_wymiar, nazwa_zbioru)
    
    return wykryty_wymiar, wariancje
end

X_hiper, _ = hiper_rozmaitosc_dane(1000)
analizuj_zbior_danych(X_hiper, "Hiper Rozmaitosc 4D w 20D", K_zakres=10:1:100, limit_wartosci=100)

X_kula, _ = kula_3d_w_15d(1000)
analizuj_zbior_danych(X_kula, "Kula 3D w 15D", K_zakres=10:1:100, limit_wartosci=100)
