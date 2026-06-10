using SparseArrays, LinearAlgebra, KrylovKit, Random
N = 1000
k = 10
W = sparse(rand(1:N, N*k), rand(1:N, N*k), rand(N*k), N, N)
M = (I - W)' * (I - W)
@time vals1, vecs1, info1 = eigsolve(M, 5, :SR; ishermitian=true, tol=1e-5)
println(vals1)
@time vals2 = sort(real(eigvals(Symmetric(Matrix(M)))))[1:5]
println(vals2)
