function update_Iρ!(qp::QP)
    c = qp.cache
    #qp.Iρ .= qp.ρ * Diagonal((qp.C * qp.x - qp.d .>= zeros(length(qp.d)) .|| qp.μ .!= 0))
    mul!(c.c_ineq, qp.C, qp.x)
    @. c.c_ineq = c.c_ineq - qp.d
    for i = 1:length(c.c_ineq)
        if (c.c_ineq[i] >= 0 || qp.μ[i] != 0)
            qp.Iρ[i, i] = qp.ρ
        else
            qp.Iρ[i, i] = 0
        end
    end
end
function update_dual!(qp::QP)
    c = qp.cache

    #qp.λ .= qp.λ .+ qp.ρ .* (qp.A * qp.x - qp.b)
    mul!(c.c_eq, qp.A, qp.x)
    @. c.c_eq = c.c_eq - qp.b
    @. qp.λ = qp.λ + qp.ρ * c.c_eq

    #qp.μ .= qp.μ + qp.Iρ * (qp.C * qp.x - qp.d)
    mul!(c.c_ineq, qp.C, qp.x)
    @. c.c_ineq = c.c_ineq - qp.d
    mul!(c.c_ineq2, qp.Iρ, c.c_ineq)
    @. qp.μ = qp.μ + c.c_ineq2
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
    mul!(c.c_ineq2, qp.Iρ, c.c_ineq)
    @. c.c_ineq2 = c.c_ineq2 + qp.μ
    mul!(c.c_n3, qp.Ct, c.c_ineq2)

    # qp.∇L .= qp.Q * qp.x + qp.q + qp.A' * (qp.λ + qp.ρ * (qp.A * qp.x - qp.b)) + qp.C' * (qp.μ + qp.Iρ * (qp.C * qp.x - qp.d))
    @. qp.∇L = c.c_n1 + qp.q + c.c_n2 + c.c_n3

    #qp.∇2L
    # qp.C' * qp.Iρ * qp.C
    mul!(c.c_ineq_n, qp.Iρ, qp.Cd)
    mul!(c.c_n_n, qp.Ct, c.c_ineq_n)

    #qp.∇2L .= qp.Q + qp.ρ .* qp.A' * qp.A + qp.C' * qp.Iρ * qp.C
    @. qp.∇2L = qp.Q + qp.ρ * qp.AtA + c.c_n_n
end
function update_penalty!(qp::QP)
    qp.ρ = qp.ρ * qp.ϕ
end
function logging(qp::QP)
    @printf("%3d   %10.3e  %9.2e  %9.2e  %9.2e\n",
        qp.iter, qp.cost, qp.eq_res, qp.ineq_res, qp.complementarity_res)
end
function check_convergence!(qp::QP)
    c = qp.cache
    #J = 0.5 * qp.x' * qp.Q * qp.x + dot(qp.q, qp.x)
    mul!(c.c_n1, qp.Q, qp.x)
    qp.cost = 0.5 * dot(qp.x, c.c_n1) + dot(qp.x, qp.q)

    #eq_res = norm(qp.A * qp.x - qp.b)
    mul!(c.c_eq, qp.A, qp.x)
    @. c.c_eq = c.c_eq - qp.b
    qp.eq_res = norm(c.c_eq)

    mul!(c.c_ineq, qp.C, qp.x)
    @. c.c_ineq = c.c_ineq - qp.d

    for i = 1:length(qp.d)
        if (c.c_ineq[i] < 0)
            c.c_ineq2[i] = -c.c_ineq[i] * qp.μ[i]
        else
            c.c_ineq2[i] = -c.c_ineq[i] * qp.μ[i]
        end
    end
    qp.complementarity_res = norm(c.c_ineq2)

    for i = 1:length(qp.d)
        if (c.c_ineq[i] < 0)
            c.c_ineq[i] = 0
        end
    end
    qp.ineq_res = norm(c.c_ineq)
    qp.converged = (qp.eq_res <= length(qp.b) * qp.tol.eq_feas) && (qp.ineq_res <= length(qp.d) * qp.tol.ineq_feas) && (qp.complementarity_res <= length(qp.d) * qp.tol.complementarity)
end