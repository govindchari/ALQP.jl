function initialize!(qp::QP)
    update_derivatives!(qp)
end
function solve!(qp::QP)
    solve!(qp, false)
end
function solve!(qp::QP, verbose::Bool)
    initialize!(qp)
    if verbose
        println("iter     objv       |Ax-b|    |Cx-d|  Complementarity\n")
        println("-------------------------------------------------------------\n")
    end
    while (!qp.converged && qp.iter < qp.tol.max_iter)
        minimize_augmented_lagrangian!(qp)
        update_dual!(qp)
        update_penalty!(qp)
        check_convergence!(qp)
        qp.iter = qp.iter + 1
        if verbose
            logging(qp)
        end
    end
    @assert qp.iter < qp.tol.max_iter
end
