using ALQP
using OSQP
using BenchmarkTools

include("portfolio.jl")

N = 10
n_list = [2^i for i = 1:N]
osqp_time = zeros(N)
alqp_time = zeros(N)
tol = 1e-3

for i = 1:N
    Q, q, A, b, C, d, A_osqp, l, u = porfolio(n_list[i])
    qp = QP(Q, q, A, b, C, d)
    m = OSQP.Model()
    OSQP.setup!(m; P=Q, q=q, A=A_osqp, l=l, u=u, verbose=false)
    result = OSQP.solve!(m)
    alqp_time[i] = @elapsed solve!(qp)
    osqp_time[i] = @elapsed OSQP.solve!(m)
    @assert norm(qp.x - result.x) < sqrt(n_list[i])*tol
end

println(n_list)
println(osqp_time)
println(alqp_time)
