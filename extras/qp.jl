using LinearAlgebra, SparseArrays
using BenchmarkTools
using Printf
using OSQP

include("portfolio.jl")

mutable struct QP
    # problem data
    Q::SparseMatrixCSC{Float64,Int64}
    q::Array{Float64,1}
    A::SparseMatrixCSC{Float64,Int64}
    b::Array{Float64,1}
    C::SparseMatrixCSC{Float64,Int64}
    d::Array{Float64,1}

    ∇L::Array{Float64,1}
    ∇2L::Array{Float64,2}
    Iρ::SparseMatrixCSC{Float64,Int64}

    x::Array{Float64,1}
    λ::Array{Float64,1}
    μ::Array{Float64,1}
    ρ::Float64
    ϕ::Float64

    aug_lag_tol::Float64

    function QP(Q, q, A, b, C, d)
        x = zeros(length(q))
        λ = zeros(length(b))
        μ = zeros(length(d))
        ρ = 1
        Iρ = I(length(d))
        ∇L = zeros(length(q))
        ∇2L = zeros(length(q), length(q))
        ϕ = 10.0
        aug_lag_tol = 1e-8
        new(Q, q, A, b, C, d, ∇L, ∇2L, Iρ, x, λ, μ, ρ, ϕ, aug_lag_tol)
    end
end
function update_dual!(qp::QP)
    qp.λ .= qp.λ .+ qp.ρ .* (qp.A * qp.x - qp.b)
    #qp.μ .= max.(0, qp.μ .+ qp.ρ .* (qp.C * qp.x - qp.d))
    qp.μ .= qp.μ + qp.Iρ * (qp.C * qp.x - qp.d)
end
function update_penalty!(qp::QP)
    qp.ρ = qp.ρ * qp.ϕ
end
function newton_step!(qp::QP)
    qp.x .= qp.x .- qp.∇2L \ qp.∇L
end
function update_Iρ!(qp::QP)
    qp.Iρ .= qp.ρ * Diagonal((qp.C * qp.x - qp.d .>= zeros(length(qp.d)) .|| qp.μ .!= 0))
end
function update_derivatives!(qp::QP)
    qp.∇L .= qp.Q * qp.x + qp.q + qp.A' * (qp.λ + qp.ρ * (qp.A * qp.x - qp.b)) + qp.C' * (qp.μ + qp.Iρ * (qp.C * qp.x - qp.d))
    qp.∇2L .= qp.Q + qp.ρ .* qp.A' * qp.A + qp.C' * qp.Iρ * qp.C
end
function minimize_augmented_lagrangian!(qp::QP)
    while (norm(qp.∇L) > qp.aug_lag_tol)
        update_derivatives!(qp)
        update_Iρ!(qp)
        newton_step!(qp)
    end
    update_derivatives!(qp)
    update_Iρ!(qp)
end
function logging(qp::QP, iter)

    J = 0.5 * qp.x' * qp.Q * qp.x + dot(qp.q, qp.x)
    gap = dot(qp.C * qp.x - qp.d, qp.μ)
    eq_res = norm(qp.A * qp.x - qp.b)


    @printf("%3d   %10.3e  %9.2e  %9.2e\n",
        iter, J, gap, eq_res)

    return nothing
end
function initialize!(qp::QP)
    update_derivatives!(qp)
    update_Iρ!(qp)
end

function solve!(qp::QP)
    solve!(qp, false)
end
function solve!(qp::QP, verbose::Bool)
    if verbose
        println("iter     objv        gap       |Ax-b|    |Gx+s-h|\n")
        println("-------------------------------------------------\n")
    end
    for i = 1:70
        minimize_augmented_lagrangian!(qp)
        update_dual!(qp)
        if (i%10 == 0)
            update_penalty!(qp)
        end
        if verbose
            logging(qp, i)
        end
        #check_convergence!(qp)
    end
end

let
    n = 1024
    tol = 1e-3
    Q, q, A, b, C, d, A_osqp, l, u = porfolio(n)
    qp = QP(Q, q[:], A, b[:], C, d[:])
    m = OSQP.Model()
    OSQP.setup!(m; P=Q, q=q, A=A_osqp, l=l, u=u, verbose=false)
    result = OSQP.solve!(m)
    @time solve!(qp,true)
    @assert norm(qp.x - result.x) < sqrt(n)*tol 
end