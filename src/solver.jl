function solve!(qp::QP)
    solve!(qp, false)
end
function solve!(qp::QP, verbose::Bool)
    if verbose
        println("iter     objv       |Ax-b|    ineq res  comp slackness\n")
        println("-------------------------------------------------------------\n")
    end
    while (!qp.converged && qp.iter < qp.tol.max_iter)
        minimize_augmented_lagrangian!(qp)
        update_dual!(qp)  
        if (qp.iter % 10 == 0)
            update_penalty!(qp)
        end
        check_convergence!(qp)
        qp.iter = qp.iter + 1
        if verbose
            logging(qp)
        end
    end
    @assert qp.iter < qp.tol.max_iter "Max Iterations Reached"
end
