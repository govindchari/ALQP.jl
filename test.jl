using ALQP
using LinearAlgebra
using SparseArrays
using OSQP
using BenchmarkTools

# ==============================KNOWN TESTCASE========================================

n = 2
Q = [1 0;0 1]
q = zeros(n)

A = zeros(n,n)
b = zeros(n)
C = [1 1]
d = [-5]

Q = sparse(Q)
A = sparse(A)
C = sparse(C)

qp = QP(Q,q,A,b,C,d)
solve!(qp, true)
println(qp.x)

# ==============================ALLOCATION TESTCASE====================================


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
m = OSQP.Model()
OSQP.setup!(m; P=Q, q=q, A=sparse(I(10)), l=zeros(n), u=ones(n), verbose=false)

#Benchmarking
println("OSQP Benchmark:")
@btime OSQP.solve!($m)

println("ALQP Benchmark:")
@btime solve!($qp)
