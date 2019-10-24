function stepOcean_Flow!(
    ocn  :: Ocean;
    cfgs...
)
    
    adv_scheme = cfgs[:adv_scheme]
    do_convadjust = cfgs[:do_convadjust]
    Δt         = cfgs[:Δt]

    if adv_scheme == :static
        return
    end

    # Determine the temperature / salinity of FLDO layer
    @loop_hor ocn i j let

        ocn.ΔT[i, j] = mixFLDO!(
            qs   = ocn.cols.Ts[i, j],
            zs   = ocn.cols.zs[i, j],
            hs   = ocn.cols.hs[i, j],
            q_ML = ocn.T_ML[i, j],
            FLDO = ocn.FLDO[i, j],
            FLDO_ratio_top = ocn.FLDO_ratio_top[i, j],
            FLDO_ratio_bot = ocn.FLDO_ratio_bot[i, j],
        )

        ocn.ΔS[i, j] = mixFLDO!(
            qs   = ocn.cols.Ss[i, j],
            zs   = ocn.cols.zs[i, j],
            hs   = ocn.cols.hs[i, j],
            q_ML = ocn.S_ML[i, j],
            FLDO = ocn.FLDO[i, j],
            FLDO_ratio_top = ocn.FLDO_ratio_top[i, j],
            FLDO_ratio_bot = ocn.FLDO_ratio_bot[i, j],
        )

    end

    # Pseudo code
    # 1. assign velocity field
    # 2. calculate temperature & salinity flux
    # 3. calculate temperature & salinity flux divergence
    # Gov eqn adv + diff: ∂T/∂t = - 1 / (ρ H1) ( ∇⋅(M1 T1) - (∇⋅M1) Tmid )
  
    calDiffAdv_QUICK!(
        ocn,
        qs          = ocn.Ts,
        wq_bnd      = ocn.wT,
        dΔqdt       = ocn.dΔTdt,
        FLUX_CONV   = ocn.TFLUX_CONV,
        FLUX_CONV_h = ocn.TFLUX_CONV_h,
        FLUX_DEN_x  = ocn.TFLUX_DEN_x,
        FLUX_DEN_y  = ocn.TFLUX_DEN_y,
        FLUX_DEN_z  = ocn.TFLUX_DEN_z,
    )

    calDiffAdv_QUICK!(
        ocn,
        qs          = ocn.Ss,
        wq_bnd      = ocn.wS,
        dΔqdt       = ocn.dΔSdt,
        FLUX_CONV   = ocn.SFLUX_CONV,
        FLUX_CONV_h = ocn.SFLUX_CONV_h,
        FLUX_DEN_x  = ocn.SFLUX_DEN_x,
        FLUX_DEN_y  = ocn.SFLUX_DEN_y,
        FLUX_DEN_z  = ocn.SFLUX_DEN_z,
    )


    @loop_hor ocn i j let
 
        Nz = ocn.Nz[i, j]
        zs   = ocn.cols.zs[i, j]
        hs   = ocn.cols.hs[i, j]
        h_ML = ocn.h_ML[i, j]
        FLDO = ocn.FLDO[i, j]
        
        for k = 1:ocn.Nz[i, j]
            ocn.Ts[k, i, j] += Δt * ocn.TFLUX_CONV[k, i, j]
            ocn.Ss[k, i, j] += Δt * ocn.SFLUX_CONV[k, i, j]
        end

        # Adjust ΔT, ΔS
        ocn.ΔT[i, j] += ocn.dΔTdt[i, j] * Δt 
        ocn.ΔS[i, j] += ocn.dΔSdt[i, j] * Δt 

#=        if (i, j) == (45, 2)
           println("ocn.TFLUX_CONV_h[1:2, 45, 2]: ", ocn.TFLUX_CONV_h[1:2, i, j])
           println("ocn.TFLUX_DEN_z[1:3, 45, 2]: ", ocn.TFLUX_DEN_z[1:3, i, j])
           println("ocn.dΔTdt[45, 2]: ", ocn.dΔTdt[i, j])
           println("ocn.dΔSdt[45, 2]: ", ocn.dΔSdt[i, j])
        end
=#

        ocn.T_ML[i, j] = unmixFLDOKeepDiff!(;
            qs   = ocn.cols.Ts[i, j],
            zs   = zs,
            hs   = hs,
            h_ML = h_ML,
            FLDO = FLDO,
            Nz   = Nz,
            Δq   = ocn.ΔT[i, j],
        )

        ocn.S_ML[i, j] = unmixFLDOKeepDiff!(;
            qs   = ocn.cols.Ss[i, j],
            zs   = zs,
            hs   = hs,
            h_ML = h_ML,
            FLDO = FLDO,
            Nz   = Nz,
            Δq   = ocn.ΔS[i, j],
        )

        OC_updateB!(ocn, i, j)

        if do_convadjust
            OC_doConvectiveAdjustment!(ocn, i, j)
        end

    end

    #=
    for k = 1:ocn.Nz_bone
        DisplacedPoleCoordinate.hadv_upwind!(
            ocn.gi,
            ocn.lays.T_hadvs[k],
            ocn.lays.u[k],
            ocn.lays.v[k],
            ocn.lays.Ts[k],
            ocn.lays.mask3[k],
        )

        DisplacedPoleCoordinate.hadv_upwind!(
            ocn.gi,
            ocn.lays.S_hadvs[k],
            ocn.lays.u[k],
            ocn.lays.v[k],
            ocn.lays.Ss[k],
            ocn.lays.mask3[k],
        )

    end

    @loop_hor ocn i j let
        
        Nz = ocn.Nz[i, j]
        zs   = ocn.cols.zs[i, j]
        hs   = ocn.cols.hs[i, j]
        h_ML = ocn.h_ML[i, j]
        FLDO = ocn.FLDO[i, j]
        Nz   = ocn.Nz[i, j]

        vadv_upwind!(
            ocn.cols.T_vadvs[i, j],
            ocn.cols.w[i, j],
            ocn.cols.Ts[i, j],
            ocn.cols.hs[i, j],
            Nz,
        )

        vadv_upwind!(
            ocn.cols.S_vadvs[i, j],
            ocn.cols.w[i, j],
            ocn.cols.Ss[i, j],
            ocn.cols.hs[i, j],
            Nz,
        )



        # Here I choose not to update the bottom layer if there is only one layer.
        ocn.wT[i, j] = ocn.cols.w[i, j][Nz] * ocn.cols.Ts[i, j][Nz]
        if (i, j) == (20, 2)
            println(format("20,2: T_ML = {:f}, T = {:f}", ocn.T_ML[20, 2], ocn.Ts[1, 20, 2] ))
            println(format("20,2: vadv = {:e}, hadv = {:e}", ocn.T_vadvs[1, 20, 2], ocn.T_hadvs[1, 20, 2] ))
        end
        for k = 1:Nz
            ocn.Ts[k, i, j] += Δt * ( ocn.T_vadvs[k, i, j] + ocn.T_hadvs[k, i, j] )
            ocn.Ss[k, i, j] += Δt * ( ocn.S_vadvs[k, i, j] + ocn.S_hadvs[k, i, j] )
            #ocn.Ts[k, i, j] += Δt * ( ocn.T_hadvs[k, i, j] )
            #ocn.Ss[k, i, j] += Δt * ( ocn.S_hadvs[k, i, j] )

        end

        #=
        if (i, j) == (80, 50)
            println("#BEFORE")
            println(ocn.cols.Ts[i,j][1:10])
            println(ocn.T_ML[i,j])
        end

        =#

        # Remix top layers
        #=
        ocn.T_ML[i, j] = remixML!(;
            qs   = ocn.cols.Ts[i, j],
            zs   = zs,
            hs   = hs,
            h_ML = h_ML,
            FLDO = FLDO,
            Nz   = Nz,
        )

        ocn.S_ML[i, j] = remixML!(;
            qs   = ocn.cols.Ss[i, j],
            zs   = zs,
            hs   = hs,
            h_ML = h_ML,
            FLDO = FLDO,
            Nz   = Nz,
        )
        =#

        if (i, j) == (20, 2)
            println(format("Before 20,2: T_ML = {:f}, T = {:f}", ocn.T_ML[20, 2], ocn.Ts[1, 20, 2] ))
        end

        ocn.T_ML[i, j] = unmixFLDOKeepDiff!(;
            qs   = ocn.cols.Ts[i, j],
            zs   = zs,
            hs   = hs,
            h_ML = h_ML,
            FLDO = FLDO,
            Nz   = Nz,
            Δq   = ocn.ΔT[i, j],
            verbose = (i, j) == (20, 2)
        )
 
        ocn.S_ML[i, j] = unmixFLDOKeepDiff!(;
            qs   = ocn.cols.Ss[i, j],
            zs   = zs,
            hs   = hs,
            h_ML = h_ML,
            FLDO = FLDO,
            Nz   = Nz,
            Δq   = ocn.ΔS[i, j],
        )
 
        if (i, j) == (20, 2)
            println(format("20,2: T_ML = {:f}, T = {:f}", ocn.T_ML[20, 2], ocn.Ts[1, 20, 2] ))
        end
      
        #=
        if (i, j) == (80, 50)
            println("#AFTER")
            println(ocn.cols.Ts[i,j][1:10])
            println(ocn.T_ML[i,j])
        end
        =#
        OC_updateB!(ocn, i, j)


        if do_convadjust
            OC_doConvectiveAdjustment!(ocn, i, j)
        end


    end
    =#

end


#=
function vadv_upwind!(
    vadvs  :: AbstractArray{Float64, 1},
    ws     :: AbstractArray{Float64, 1},
    qs     :: AbstractArray{Float64, 1},
    #Δzs    :: AbstractArray{Float64, 1},
    hs    :: AbstractArray{Float64, 1},
    Nz     :: Integer,
)

    #
    # Array information:
    # 
    # length(ws)     == length(qs) or Nz + 1
    # length(Δzs)    == length(qs) or Nz - 1
    # length(qstmp)  == length(qs) or Nz
    # length(flxtmp) == length(qs) or Nz
    # 
    # Nz reveals if there is bottom of ocean
    # 

    # Extreme case: only one grid point
    if Nz <= 1
        vadvs .= 0
        return
    end

#=
    if ws[1] > 0.0
        vadvs[1] = - ws[1] * (qs[1] - qs[2]) / Δzs[1]
    else
        vadvs[1] = 0.0
    end

    for k = 2:Nz-1
        if ws[k] > 0.0
            vadvs[k] = - ws[k] * (qs[k] - qs[k+1]) / Δzs[k]
        else
            vadvs[k] = - ws[k] * (qs[k-1] - qs[k]) / Δzs[k-1]
        end
    end
=#
#= 2019/10/19 commented
    if ws[1] > 0.0
        vadvs[1] = - ws[1] * (qs[1] - qs[2]) / hs[1]
    else
        vadvs[1] = 0.0
    end

    for k = 2:Nz-1
        if ws[k] > 0.0
            vadvs[k] = - ws[k] * (qs[k] - qs[k+1]) / hs[k]
        else
            vadvs[k] = - ws[k] * (qs[k-1] - qs[k]) / hs[k]
        end
    end
=#

    # Skip the last layer
#    vadvs[Nz] = 0.0
#=
    if ws[Nz] > 0.0
        vadvs[Nz] = 0.0
    else
        vadvs[Nz] = - ws[Nz] * (qs[Nz-1] - qs[Nz]) / hs[Nz]
    end
=#

    vadvs[1] = 0.0

    if ws[2] > 0.0
        vadvs[2] = - ws[2] * qs[2] / 2.0 / hs[2]

    for k = 2:Nz-1
        if ws[k] > 0.0
            vadvs[k] = - ( ((k==2) ? 0.0 : ws[k-1] * qs[k-1] ) - qs[k+1]) / hs[k]
        else
            vadvs[k] = - ws[k] * (qs[k-1] - qs[k]) / hs[k]
        end
    end

end

=#
"""

Calculation follows "MacCormack method" of Lax-Wendroff scheme.

Reference website:
  https://en.wikipedia.org/wiki/Lax%E2%80%93Wendroff_method#MacCormack_method

"""
function doZAdvection_MacCormack!(;
    Nz     :: Integer,
    qs     :: AbstractArray{Float64, 1},
    ws     :: AbstractArray{Float64, 1},
    Δzs    :: AbstractArray{Float64, 1},
    Δt     :: Float64,
    qstmp  :: AbstractArray{Float64, 1},
    flxtmp :: AbstractArray{Float64, 1},
)

    #
    # Array information:
    # 
    # length(ws)     == length(qs) or Nz + 1
    # length(Δzs)    == length(qs) or Nz - 1
    # length(qstmp)  == length(qs) or Nz
    # length(flxtmp) == length(qs) or Nz
    # 
    # Nz reveals if there is bottom of ocean
    # 

    # Step 1: calculate flux
    for k = 1:Nz
        flxtmp[k] = (ws[k] + ws[k+1]) / 2.0  * qs[k]
    end
   
    # Step 2: Calculate qs using flux from bottom
    for k = 1:length(qstmp)-1
        qstmp[k] = qs[k] - Δt / Δzs[k] * (flxtmp[k] - flxtmp[k+1])
    end
    # ignore flux at the last layer

    # Step 3: calculate updated flux
    for k = 1:Nz
        flxtmp[k] = (ws[k] + ws[k+1]) / 2.0  * qstmp[k]
    end
 
    # Step 4: Flux from top
    
    # First layer has no flux from above
    qs[1] = (qstmp[1] + qs[1]) / 2.0  - Δt / Δzs[1] / 2.0 * (- flxtmp[1])
    for k = 2:length(qstmp)-1
        qs[k] = (qstmp[k] + qs[k]) / 2.0 - Δt / Δzs[k-1] / 2.0 * (flxtmp[k-1] - flxtmp[k])
    end
   
    # Don't update qs[end] 
    
end

