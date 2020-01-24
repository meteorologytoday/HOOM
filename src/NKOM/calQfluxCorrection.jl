
#=
function doFluxCorrection!(;
    qs        :: AbstractArray{Float64, 2},
    qs_target :: AbstractArray{Float64, 2},
    Δt        :: Float64,
    τ         :: Float64 = 15 * 86400,
)



end
=#

function calQflux_correction!(
    ocn :: Ocean;
    τ   :: Float64 = 15 * 86400.0,
    cfgs...
)

    Δt = cfgs[:Δt]
    r = Δt / τ
    rr = r / (1.0 + r)

    @loop_hor ocn i j let
       
        # Euler backward method

        if (i, j) == (50, 50)
            println("ocn.in_flds.sst[i, j] = ", ocn.in_flds.sst[i, j])
            println("ocn.T_ML[i, j] = ", ocn.T_ML[i, j])
            println("ΔT = ", r * (ocn.in_flds.sst[i, j] - ocn.T_ML[i, j]) / (1.0 + r))
        end


        T_ML = ocn.T_ML[i, j]
        FLDO = ocn.FLDO[i, j]
 
        ΔT = rr * (ocn.in_flds.sst[i, j] - T_ML)
        T_ML += ΔT
        ocn.T_ML[i, j] = T_ML
        if FLDO > 1
            ocn.Ts[1:FLDO-1, i, j] .= T_ML
        elseif FLDO == -1
            ocn.Ts[1:ocn.Nz[i, j], i, j] .= T_ML
        end

        ocn.qflx_correction[i, j] = - ΔT * ocn.h_ML[i, j] * ρc / Δt   # neg => warming

        if (i, j) == (50, 50)
            println("qflx_correction = ",  ocn.qflx_correction[i, j])
        end


    end

end
