
#=
function doFluxCorrection!(;
    qs        :: AbstractArray{Float64, 2},
    qs_target :: AbstractArray{Float64, 2},
    Δt        :: Float64,
    τ         :: Float64 = 15 * 86400,
)



end
=#

function calFlxCorrection!(
    ocn :: Ocean;
    τ   :: Float64 = 15 * 86400.0,
    cfgs...
)
    do_convadjust = cfgs[:do_convadjust]

    Δt = cfgs[:Δt]
    r = Δt / τ
    rr = r / (1.0 + r)



    @loop_hor ocn i j let
       
        # Euler backward method

        ifrac = ocn.in_flds.ifrac[i, j]
        ifrac_clim = ocn.in_flds.IFRACclim[i, j]
        T_ML = ocn.T_ML[i, j]
        S_ML = ocn.S_ML[i, j]
        FLDO = ocn.FLDO[i, j]
        h_ML = ocn.h_ML[i, j]
 

        ΔS = rr * (ocn.in_flds.Sclim[i, j] - S_ML)

        # when ifrac == 0 then SST is not changing while energy are still flowing
        # the solution is to use ifrac as additional information to recover Qflx_T
        ΔT_openocn = rr * (ocn.in_flds.Tclim[i, j] - T_ML)

        # Old method. Keep in comment for safety
        # Assume seaice thickenss = 1.0m
        # energy_to_melt_seaice = 1.0 * (ifrac - ifrac_clim) * ρ_si * Hf_sw
        # ΔT_seaice = 100energy_to_melt_seaice / h_ML / ρc_sw * r

        # New method: observe that between -1~-2 degC, slope of SST - IFRAC is
        # roughly 100% / 1K. Use this as a diagnostic relation to determine
        # nudged SST.
        ΔT_seaice = rr * (ifrac - ifrac_clim)

        # If ifrac is too high, then the melting will be intensified
        # Reason is that there is an upper bound for ifrac that limits
        # the relaxation rate.
        if ifrac > .90 && ΔT_seaice > 0.0
            ΔT_seaice *= 1.0 + 9.0 * (ifrac - .90) / (1.00 - .90)
        end

        #=
        IFRAC
         0.30 => 1  => dominated by ΔT_seaice
         0.15 => 0  => dominated by ΔT_openocn
        =#
        wgt_seaice= max( min( (1 - 0) / (.30 - .15) * (ifrac - .15), 1.0), 0.0)
        ΔT = ΔT_openocn * (1.0 - wgt_seaice) + ΔT_seaice * wgt_seaice 


        T_ML += ΔT
        S_ML += ΔS
        ocn.T_ML[i, j] = T_ML
        ocn.S_ML[i, j] = S_ML
        if FLDO > 1
            ocn.Ts[1:FLDO-1, i, j] .= T_ML
            ocn.Ss[1:FLDO-1, i, j] .= S_ML
        elseif FLDO == -1
            ocn.Ts[1:ocn.Nz[i, j], i, j] .= T_ML
            ocn.Ss[1:ocn.Nz[i, j], i, j] .= S_ML
        end

        ocn.qflx_T_correction[i, j] = ΔT * ocn.h_ML[i, j] * ρc_sw   / Δt   # + => warming
        ocn.qflx_S_correction[i, j] = ΔS * ocn.h_ML[i, j]           / Δt   # + => saltier

        OC_updateB!(ocn, i, j)

    end

    if do_convadjust
        @loop_hor ocn i j let
            OC_doConvectiveAdjustment!(ocn, i, j;)
        end
    end


end
