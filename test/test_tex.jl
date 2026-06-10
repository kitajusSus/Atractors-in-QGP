using LaTeXStrings
k_values = 5:5:50
k_min = minimum(k_values)
k_max = maximum(k_values)
τ = 0.25
str = L"\text{LLE analiza wartości własnych (dla K z zakresu (%$(k_min) - %$(k_max)) (średnia } \pm \sigma \text{) dla } \tau=%$τ"
println(str)
