using LinearAlgebra
using SparseArrays

function porfolio(n::Int64)
    Q = rand(n, n)
    Q = Q' * Q + n * I(n)
    Q = sparse(Q)
    q = -rand(n, 1)
    A = ones(1, n)
    b = [1.0]
    C = -I(n)
    d = zeros(n, 1)
    A_osqp = zeros(n + 1, n)
    A_osqp[1:n, 1:n] .= I(n)
    A_osqp[n+1, :] = ones(1, n)
    A_osqp = sparse(A_osqp)
    l = zeros(n + 1, 1)
    l[n+1] = 1.0
    u = ones(n + 1, 1)
    u = typemax(Float64) * u
    u[n+1] = 1.0
    return Q, q[:], A, b, C, d, A_osqp, l[:], u[:]    
end