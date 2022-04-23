function newton_step!(qp::QP)
    qp.x .= qp.x .- qp.∇2L \ qp.∇L
end
function update_derivatives!(qp::QP)
    qp.∇L .= qp.Q * qp.x + qp.q + qp.A' * (qp.λ + qp.ρ * (qp.A * qp.x - qp.b)) + qp.C' * (qp.μ + qp.Iρ * (qp.C * qp.x - qp.d))
    qp.∇2L .= qp.Q + qp.ρ .* qp.A' * qp.A + qp.C' * qp.Iρ * qp.C
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