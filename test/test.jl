using ALQP
using LinearAlgebra
using SparseArrays
using OSQP
using BenchmarkTools

include("portfolio.jl")

tol = 2e-4

# ==============================KNOWN TESTCASE========================================
n = 2
Q = [1 0; 0 1]
q = zeros(n)

A = zeros(n, n)
b = zeros(n)
C = [1 1]
d = [-5]

Q = sparse(Q)
A = sparse(A)
C = sparse(C)

qp = QP(Q, q, A, b, C, d)
solve!(qp)
@assert norm(qp.x - [-2.5; -2.5]) < sqrt(n)*tol

# ==============================PORTFOLIO TESTCASE====================================
n = 1000
Q, q, A, b, C, d, A_osqp, l, u = porfolio(n)
qp = QP(Q, q, A, b, C, d)
m = OSQP.Model()
OSQP.setup!(m; P=Q, q=q, A=A_osqp, l=l, u=u, verbose=false)
result = OSQP.solve!(m)
solve!(qp)
@assert norm(qp.x - result.x) < sqrt(n)*tol

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
solve!(qp)
result = OSQP.solve!(m)
@assert norm(qp.x - result.x) < sqrt(n)*tol
