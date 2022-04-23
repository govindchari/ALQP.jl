using ALQP
using LinearAlgebra
using SparseArrays

n = 10
Q = randn(n, n)
Q = Q' * Q
Q = sparse(Q)
q = zeros(n)
A = spzeros(0, n)
b = []
G = sparse([I(n); -I(n)])
h = [ones(n); zeros(n)]

qp = QP(Q, q, A, b, G, h)

solve!(qp)  


#=
m = OSQP.Model()
OSQP.setup!(m; P=Q, q=q, A=sparse(I(10)), l=zeros(n), u=ones(n), verbose=false)

#Benchmarking
println("OSQP Benchmark:")
@btime OSQP.solve!($m)

println("ALQP Benchmark:")
@btime solve!($qp)
=#
