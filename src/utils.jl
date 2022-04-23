function update_Iρ!(qp::QP)
    qp.Iρ .= qp.ρ * Diagonal((qp.C * qp.x - qp.d .>= zeros(length(qp.d)) .|| qp.μ .!= 0))
end
function update_derivatives!(qp::QP)
    qp.∇L .= qp.Q * qp.x + qp.q + qp.A' * (qp.λ + qp.ρ * (qp.A * qp.x - qp.b)) + qp.C' * (qp.μ + qp.Iρ * (qp.C * qp.x - qp.d))
    qp.∇2L .= qp.Q + qp.ρ .* qp.A' * qp.A + qp.C' * qp.Iρ * qp.C
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