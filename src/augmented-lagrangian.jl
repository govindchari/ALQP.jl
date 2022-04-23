function newton_step!(qp::QP)
    qp.cache.c_n1 .= qp.∇2L \ qp.∇L
    @. qp.x = qp.x - qp.cache.c_n1
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
    mul!(c.c_n2, At, c.c_eq)

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
    mul!(c.c_n_n, Ct, qp.Iρ)
    mul!(c.c_n_n, c.c_n_n, qp.C)

    #qp.∇2L .= qp.Q + qp.ρ .* qp.A' * qp.A + qp.C' * qp.Iρ * qp.C
    @. qp.∇2L = qp.Q + qp.ρ * qp.AtA + c.c_n_n
end
function minimize_augmented_lagrangian!(qp::QP)
    while (norm(qp.∇L) > qp.tol.aug_lag)
        update_derivatives!(qp)
        update_Iρ!(qp)
        newton_step!(qp)
    end
    update_derivatives!(qp)
    update_Iρ!(qp)
end