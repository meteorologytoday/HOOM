function stepColumn!(
    mb :: ModelBlock,
    Δt :: Float64,
)

    fi = mb.fi
    co = mb.co

    # *** Diffusion and restoring ***
    # (b_t+1 - b_t) / dt =  OP1 * b_t+1 + OP2 * (b_t+1 - b_target) + const
    # b_t+1 - b_t =  dt * OP1 * b_t+1 + dt * OP2 * (b_t+1 - b_target) + dt * const
    # (I - OP1 * dt - OP2 * dt) b_t+1 = b_t - dt * OP2 * b_target + dt * const
    # b_t+1 = (I - OP1 * dt - OP2 * dt) \ (b_t - dt * OP2 * b_target + dt * const)

    op_vdiff = calOp_vdiff(co.vd, view(fi._X_, :, 1), view(fi.HMXL, :))

    op = op_vdiff
    F_EBM = lu( I - Δt * op)
    rad = ( co.mtx[:T_swflxConv_sT] * view(fi.SWFLX, :) + co.mtx[:T_nswflxConv_sT] * view(fi.NSWFLX, :)) / ρcp_sw
    #fi._X_[:, 1] = F_EBM \ ( fi._X_[:, 1] + Δt * rad)
    fi._X_[:, 1] =  ( fi._X_[:, 1] + Δt * rad)
    
end
