
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
        T_ML = ocn.T_ML[i, j]
        S_ML = ocn.S_ML[i, j]
        FLDO = ocn.FLDO[i, j]
        h_ML = ocn.h_ML[i, j]
 
        ΔT = rr * (ocn.in_flds.Tclim[i, j] - T_ML)
        ΔS = rr * (ocn.in_flds.Sclim[i, j] - S_ML)

        # Assume seaice thickenss = 1.0m
        energy_to_melt_seaice = 1.0 * (ifrac - ocn.in_flds.IFRACclim[i, j]) * ρ_si * Hf_sw
        ΔT_to_melt_seaice = energy_to_melt_seaice / h_ML / ρc_sw * r

        ΔT += ΔT_to_melt_seaice

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
