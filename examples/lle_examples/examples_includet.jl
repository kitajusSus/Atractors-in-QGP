include("helix.jl")
include("swissroll.jl")
include("scurve.jl")
include("2d_examples_data.jl")
include("../local_pca/lpca.jl")


function symulacja_dims(N,K)


    zbiory_danych = [
    ("Płaski Dywan", flat_carpet_dane),
    ("Swiss Roll", swissroll_dane),
    ("Krzywa S", scurve_dane),
    ("Sinus", sinus_dane),
    ("Hiperbola", hiperbola_dane)
    ]

    for (nazwa, funkcja) in zbiory_danych
        X, _ = funkcja(N)
        W = ncbj4_lle_basic(X, K)
        M = (I - W)' * (I - W)
        F = eigen(Symmetric(M))
        
        println(nazwa)
        println("----------")
        println(F.values)
    
    end

    for (nazwa, funkcja) in zbiory_danych
        X, _ = funkcja(N)
        dd = dims(X, k=K)
        meandim = mean(dd)
        println("$nazwa - dla K = $(K) Średnia wymiarowość: $meandim")
    end

end
