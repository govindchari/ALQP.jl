function initialize!(qp::QP)
    update_derivatives!(qp)
end
function solve!(qp::QP)
    solve!(qp, false)
end
function solve!(qp::QP, verbose::Bool)
    initialize!(qp)
    if verbose
        println("iter     objv       |Ax-b|    ineq res  comp slackness\n")
        println("-------------------------------------------------------------\n")
    end
    while (!qp.converged && qp.iter < qp.tol.max_iter)
        if (qp.iter<=80 && qp.iter%10 == 0)
            qp.ρ = 10^(qp.iter/10)
        end
        minimize_augmented_lagrangian!(qp)
        update_dual!(qp)
        update_penalty!(qp)
        check_convergence!(qp)
        qp.iter = qp.iter + 1
        if verbose
            logging(qp)
        end
    end
    @assert qp.iter < qp.tol.max_iter "Max Iterations Reached"
end
