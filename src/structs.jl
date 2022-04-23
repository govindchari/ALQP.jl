struct TOLERANCE
    aug_lag::Float64
    eq_feas::Float64
    ineq_feas::Float64
    complementarity::Float64
    max_iter::Int64

    function TOLERANCE()
        aug_lag = 1e-8
        eq_feas = 1e-6
        ineq_feas = 1e-6
        complementarity = 1e-2
        max_iter = 100
        new(aug_lag, eq_feas, ineq_feas, complementarity, max_iter)
    end
end

struct QP
    Q::SparseMatrixCSC{Float64,Int64}
    q::Array{Float64,1}
    A::SparseMatrixCSC{Float64,Int64}
    b::Array{Float64,1}
    C::SparseMatrixCSC{Float64,Int64}
    d::Array{Float64,1}

    At::SparseMatrixCSC{Float64,Int64}
    Ct::SparseMatrixCSC{Float64,Int64}

    ∇L::Array{Float64,1}
    ∇2L::Array{Float64,2}
    Iρ::SparseMatrixCSC{Float64,Int64}
    AtA::SparseMatrixCSC{Float64,Int64}

    x::Array{Float64,1}
    λ::Array{Float64,1}
    μ::Array{Float64,1}
    ρ::Float64
    ϕ::Float64

    tol::TOLERANCE
    cache::CACHE

    function QP(Q, q, A, b, C, d)
        n = length(q)
        neq = length(b)
        nineq = length(d)
        x = zeros(length(q))
        λ = zeros(length(b))
        μ = zeros(length(d))
        ρ = 1
        Iρ = I(length(d))
        AtA = A' * A
        ∇L = zeros(length(q))
        ∇2L = zeros(length(q), length(q))
        ϕ = 1.0

        tol = TOLERANCE()
        cache = CACHE(neq, nineq, n)
        new(Q, q, A, b, C, d, A', C', ∇L, ∇2L, Iρ, AtA, x, λ, μ, ρ, ϕ, tol, cache)
    end
end
