function update_Iρ!(qp::QP)
    qp.Iρ .= qp.ρ * Diagonal((qp.C * qp.x - qp.d .>= zeros(length(qp.d)) .|| qp.μ .!= 0))
end
function update_dual!(qp::QP)
    #qp.λ .= qp.λ .+ qp.ρ .* (qp.A * qp.x - qp.b)
    qp.μ .= max.(0, qp.μ .+ qp.ρ .* (qp.C * qp.x - qp.d))
    #qp.μ .= qp.μ + qp.Iρ * (qp.C * qp.x - qp.d)
end
function update_derivatives!(qp::QP)
    c = qp.cache

    # qp.∇L
    # Q*x
    mul!(c.c_n1, qp.Q, qp.x)

    # qp.A' * (qp.λ + qp.ρ * (qp.A * qp.x - qp.b))
    mul!(c.c_eq, qp.A, qp.x)
    @. c.c_eq = c.c_eq - qp.b
    @. c.c_eq = qp.ρ * c.c_eq
    @. c.c_eq = c.c_eq + qp.λ
    mul!(c.c_n2, qp.At, c.c_eq)

    # qp.C' * (qp.μ + qp.Iρ * (qp.C * qp.x - qp.d))
    mul!(c.c_ineq, qp.C, qp.x)
    @. c.c_ineq = c.c_ineq - qp.d
    mul!(c.c_ineq, qp.Iρ, c.c_ineq)
    @. c.c_ineq = c.c_ineq + qp.μ
    mul!(c.c_n3, qp.Ct, c.c_ineq)

    # qp.∇L .= qp.Q * qp.x + qp.q + qp.A' * (qp.λ + qp.ρ * (qp.A * qp.x - qp.b)) + qp.C' * (qp.μ + qp.Iρ * (qp.C * qp.x - qp.d))
    @. qp.∇L = c.c_n1 + c.c_n2 + c.c_n3

    #qp.∇2L
    # qp.C' * qp.Iρ * qp.C
    mul!(c.c_ineq_n, qp.Iρ, qp.C)
    mul!(c.c_n_n, qp.Ct, c.c_ineq_n)

    #qp.∇2L .= qp.Q + qp.ρ .* qp.A' * qp.A + qp.C' * qp.Iρ * qp.C
    @. qp.∇2L = qp.Q + qp.ρ * qp.AtA + c.c_n_n
end
function update_penalty!(qp::QP)
    qp.ρ = qp.ρ * qp.ϕ
end
function logging(qp::QP, iter)

    J = 0.5 * qp.x' * qp.Q * qp.x + dot(qp.q, qp.x)
    gap = dot(qp.C * qp.x - qp.d, qp.μ)
    eq_res = norm(qp.A * qp.x - qp.b)


    @printf("%3d   %10.3e  %9.2e  %9.2e\n",
        iter, J, gap, eq_res)

end