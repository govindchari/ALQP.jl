function newton_step!(qp::QP)
    qp.cache.c_n1 .= qp.∇2L \ qp.∇L
    @. qp.x = qp.x - qp.cache.c_n1
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