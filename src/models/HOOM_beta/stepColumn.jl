function stepColumn!(
    mb :: ModelBlock,
    Δt :: Float64,
)

    fi = mb.fi
    co = mb.co
    cfg = mb.ev.config

    # *** Diffusion and restoring ***
    # (b_t+1 - b_t) / dt =  OP1 * b_t+1 + OP2 * (b_t+1 - b_target) + const
    # b_t+1 - b_t =  dt * OP1 * b_t+1 + dt * OP2 * (b_t+1 - b_target) + dt * const
    # (I - OP1 * dt - OP2 * dt) b_t+1 = b_t - dt * OP2 * b_target + dt * const
    # b_t+1 = (I - OP1 * dt - OP2 * dt) \ (b_t - dt * OP2 * b_target + dt * const)

    op_vdiff = calOp_vdiff(co.vd, view(fi._X_, :, 1), fi.HMXL)

    op = op_vdiff
    
    F_EBM = lu( I - Δt * op)
    rad = ( co.mtx[:T_swflxConv_sT] * view(fi.SWFLX, :) + co.mtx[:T_nswflxConv_sT] * view(fi.NSWFLX, :)) / ρcp_sw

    RHS = fi._X_[:, 1] + Δt * rad

    if cfg[:weak_restoring] == :on
        op   += co.mtx[:T_invτwk_TEMP_T]
        RHS .-= Δt * co.mtx[:T_invτwk_TEMP_T] * view( co.datastream["TEMP"] , :)
    end
   
 
    F_EBM = lu( I - Δt * op)
    fi._X_[:, 1] = F_EBM \ RHS
    #( fi._X_[:, 1] + Δt * rad)
    #fi._X_[:, 1] = view(co.datastream["TEMP"], :)
    #fi._X_[:, 1] =  ( fi._X_[:, 1] + Δt * rad)
    
end
