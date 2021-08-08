function stepColumn!(
    mb :: ModelBlock,
    Δt :: Float64,
)

    fi = mb.fi
    tmpfi = mb.tmpfi
    co = mb.co
    cfg = mb.ev.config

    # *** Diffusion and restoring ***
    # (b_t+1 - b_t) / dt =  OP1 * b_t+1 + OP2 * (b_t+1 - b_target) + const
    # b_t+1 - b_t =  dt * OP1 * b_t+1 + dt * OP2 * (b_t+1 - b_target) + dt * const
    # (I - OP1 * dt - OP2 * dt) b_t+1 = b_t - dt * OP2 * b_target + dt * const
    # b_t+1 = (I - OP1 * dt - OP2 * dt) \ (b_t - dt * OP2 * b_target + dt * const)
    #
    # Change of b
    #
    # Δb = b_t+1 - b_t = Δt * OP1 * b_t+1 + Δt * OP2 * (b_t+1 - b_target) + Δt * const
    #
    # Notice that I carefully use the buoyancy before
    # stepAdvection which means operators should be at
    # t = t0. I empirically know if we use newton method
    # to find a steady state, then this scheme seem to
    # ensure that physical stepping will be steady state
    # too.
    # 
    op_vdiff = calOp_vdiff(co.vd, fi._b, fi.HMXL)

    op_TEMP = op_vdiff
    op_SALT = op_vdiff
  
    # Save this operator for diagnostic purpose 
    mb.tmpfi.op_vdiff = op_vdiff
 
    # Temperature
    rad = ( co.mtx[:T_swflxConv_sT] * view(fi.SWFLX, :) + co.mtx[:T_sfcflxConv_sT] * view(fi.NSWFLX, :)) / ρcp_sw
    RHS_TEMP = view(tmpfi._INTMX_, :, 1) + Δt * rad

    # Salinity
    RHS_SALT = view(tmpfi._INTMX_, :, 2) + Δt * co.mtx[:T_sfcflxConv_sT] * view(fi.VSFLX, :)

    if cfg[:weak_restoring] == :on
        op_TEMP   += co.mtx[:T_invτwk_TEMP_T]
        op_SALT   += co.mtx[:T_invτwk_SALT_T]

        RHS_TEMP .-= Δt * co.mtx[:T_invτwk_TEMP_T] * reshape( co.datastream["TEMP"] , :)
        RHS_SALT .-= Δt * co.mtx[:T_invτwk_SALT_T] * reshape( co.datastream["SALT"] , :)
    end
 
    if cfg[:Qflx] == :on
        RHS_TEMP .+= Δt * view( co.datastream["Qflx_T"] , :) / ρcp_sw
        RHS_SALT .+= Δt * view( co.datastream["Qflx_S"] , :) 
    end
   
    F_EBM_TEMP = lu( I - Δt * op_TEMP )
    F_EBM_SALT = lu( I - Δt * op_SALT )

    tmpfi.sv[:NEWTEMP][:] = F_EBM_TEMP \ RHS_TEMP
    tmpfi.sv[:NEWSALT][:] = F_EBM_SALT \ RHS_SALT
    

    # Recompute source and sink of tracers
    if cfg[:weak_restoring] == :on
        tmpfi._WKRSTΔX_[:, 1] = tmpfi._NEWX_[:, 1] - reshape(co.datastream["TEMP"], :)
        tmpfi._WKRSTΔX_[:, 2] = tmpfi._NEWX_[:, 2] - reshape(co.datastream["SALT"], :)
        fi._WKRST_[:, 1] .= co.mtx[:T_invτwk_TEMP_T] * view(tmpfi._WKRSTΔX_, :, 1)
        fi._WKRST_[:, 2] .= co.mtx[:T_invτwk_SALT_T] * view(tmpfi._WKRSTΔX_, :, 2)
    end
 
     
end
