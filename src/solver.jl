function initialize!()
    update_derivatives!(qp)
    update_Iρ!(qp)
end
function solve!(qp::QP)
    solve!(qp, false)
end
function solve!(qp::QP, verbose::Bool)
    if verbose
        println("iter     objv        gap       |Ax-b|    |Gx+s-h|\n")
        println("-------------------------------------------------\n")
    end
    for i = 1:10
        minimize_augmented_lagrangian!(qp)
        update_dual!(qp)
        update_penalty!(qp)
        if verbose
            logging(qp, i)
        end
    end
end
