function newton_step!(qp::QP)
    qp.cache.c_n1 .= qp.∇2L \ qp.∇L
    @. qp.x = qp.x - qp.cache.c_n1
end
function minimize_augmented_lagrangian!(qp::QP)
    qp.aug_lag_iter = 0
    while (norm(qp.∇L) > qp.tol.aug_lag)
        update_derivatives!(qp)
        update_Iρ!(qp)
        newton_step!(qp)
        qp.aug_lag_iter = qp.aug_lag_iter + 1
        @assert qp.aug_lag_iter < qp.tol.max_aug_lag_iter "Max Newton Iterations Reached"
    end
    update_derivatives!(qp)
    update_Iρ!(qp)
end