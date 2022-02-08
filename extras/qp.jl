using LinearAlgebra
using SparseArrays

struct QP
    # problem data
    Q::Array{Float64,2}
    q::Array{Float64,1}
    A::Array{Float64,2}
    b::Array{Float64,1}
    C::Array{Float64,2}
    d::Array{Float64,1}
    AAt::Array{Float64,2}

    x::Array{Float64,1}
    λ::Array{Float64,1}
    μ::Array{Float64,1}
    ρ::Float64

    function QP(Q,q,A,b,C,d)
        AAt = A*A'
        x = zeros(length(q))
        λ = zeros(length(b))
        μ = zeros(length(d))
        new(Q,q,A,b,C,d,AAt,x,λ,μ)
    end
end
function update_dual!(qp::QP)
    
end
function update_penalty!(qp::QP)
    
end
function newton_step!(qp::QP)

end
function linesearch!(qp::QP)
    
end
function minimize_augmented_lagrangian!(qp::QP)
    
end
function solve!(qp::QP)
    minimize_augmented_lagrangian!(qp)
    update_dual!(qp)
    update_penalty!(qp)
end

let
    
end
