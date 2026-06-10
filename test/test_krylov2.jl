using SparseArrays, LinearAlgebra, KrylovKit, Random

N = 5000
k = 10
I_idx = Int[]
J_idx = Int[]
V = Float64[]
for i in 1:N
    idx = rand(1:N, k)
    append!(I_idx, fill(i, k))
    append!(J_idx, idx)
    append!(V, rand(k))
end
W = sparse(I_idx, J_idx, V, N, N)
M = (I - W)' * (I - W)

println("Standard KrylovKit :SR (will likely fail):")
@time vals1, vecs1, info1 = eigsolve(M, 5, :SR; ishermitian=true)
println(vals1)

println("\nKrylovKit with Shift-and-Invert:")
sigma = -1e-5
F = cholesky(Symmetric(M - sigma * I))
inv_f = x -> F \ x
@time vals_inv, vecs2, info2 = eigsolve(inv_f, N, 5, :LM; ishermitian=true)
vals2 = (1 ./ vals_inv) .+ sigma
println(sort(vals2))
