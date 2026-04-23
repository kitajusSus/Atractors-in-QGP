using Random

function flat_carpet_dane(N; )
    rnd = Xoshiro(5)
    x = rand(rnd, Float64, N)
    y = rand(rnd, Float64, N)
    X = [x y]'
    
    labels = x
    
    return X, labels

end





function hiperbola_dane(N; )
    rnd = Xoshiro(5)
    x = rand(rnd, Float64, N)
    y = 1 ./x 
    hip  = [x y]'
    
    labels = x
    
    return hip, labels
end



function sinus_dane(N)
    rnd = Xoshiro(5)
    x = 10pi .* rand(rnd, Float64, N)
    y = sin.(x)
    return [x y]', x
end

function porownaj_zbiory(; N=2000, K=12, indeksy=1:6)
    zbiory = [
        ("Plaski Dywan", flat_carpet_dane),
        ("Hiperbola", hiperbola_dane),
        ("Sinus", sinus_dane)
    ]
    
    for (nazwa, funkcja) in zbiory
        X, _ = funkcja(N)
        W = ncbj4_lle_basic(X, K)
        M = (I - W)' * (I - W)
        F = eigen(Symmetric(M))
        
        limit = min(length(F.values), maximum(indeksy))
        idx = indeksy[indeksy .<= limit]
        
        wartosci = F.values[idx]
        znormalizowane = wartosci ./ wartosci[2]
        
        println(nazwa)
        println("----------")
        println(typeof(znormalizowane))
        for i in 1:length(idx)
            println(idx[i], ": ",znormalizowane[i])
        end
        println("\n")
    end
end

porownaj_zbiory(N=2000, K=12, indeksy=1:6)
