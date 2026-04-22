include("eigen_value.jl")



function generuj_tabele()
    zakres_K = 1:30
    ile_𝛌 = 2:6 
    𝐱 = zeros(Float64, length(zakres_K), length(ile_𝛌))

    for i in zakres_K
        𝐱[i, :] = ex_eigen_value(N=2000, K=i, indeksy=ile_𝛌)
    end

    return 𝐱 
end
