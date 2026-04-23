using LinearAlgebra
include("../ncbj_lle.jl")
using GLMakie
include("swissroll.jl")

function ex_eigen_value(; N=2000, K=40, indeksy = 2:10)
    X, labels = swissroll_dane(N)
    
    W = ncbj4_lle_basic(X, K)

    M = (I - W)' * (I - W)
    F = eigen(Symmetric(M))
    wartosci_wlasne = F.values
    wartosci_do_wykresu = wartosci_wlasne[indeksy]
    return wartosci_do_wykresu #./ wartosci_do_wykresu[1]
    # fig = Figure()
    # ax = Axis(fig[1, 1], title = "Wartości własne", xlabel = "Indeks", ylabel = "Wartosc")
    # scatterlines!(ax, indeksy, wartosci_do_wykresu, markersize=10)
    # display(fig)
    # return W, F.values, wartosci_do_wykresu, wartosci_wlasne

end



ex_eigen_value()
