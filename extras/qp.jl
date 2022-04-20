using LinearAlgebra

struct QP
    # problem data
    Q::Array{Float64,2}
    q::Array{Float64,1}
    A::Array{Float64,2}
    b::Array{Float64,1}
    C::Array{Float64,2}
    d::Array{Float64,1}
    AAt::Array{Float64,2}

    ∇f::Array{Float64,1}
    ∇2f::Array{Float64,2}
    α::Float64
   
    x::Array{Float64,1}
    Δx::Array{Float64,1}
    λ::Array{Float64,1}
    μ::Array{Float64,1}
    ρ::Float64
    ϕ::Float64

    function QP(Q,q,A,b,C,d)
        AAt = A*A'
        x = zeros(length(q))
        λ = zeros(length(b))
        μ = zeros(length(d))
        ∇f = zeros(length(q))
        ∇2f = zeros(length(q),length(q))
        ϕ = 10.0
        new(Q,q,A,b,C,d,AAt,∇f,∇2f,x,λ,μ,ϕ)
    end
end
function update_dual!(qp::QP)
    λ .= λ .+ ρ*(A*x-b)
    μ .= max.(0,μ+ρ*(C*x-d))
end
function update_penalty!(qp::QP)
    ρ = ρ*ϕ
end
function newton_step!(qp::QP)
    Δx .= -∇2f\∇f
    x .= x.+Δx
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
