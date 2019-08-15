using Statistics

function stepOcean_prepare!(ocn::Ocean; cfgs...)

    adv_scheme = cfgs[:adv_scheme]

    if adv_scheme == :static
        return
    end

    # Transform input wind stress vector first
    DisplacedPoleCoordinate.project!(ocn.gi, ocn.in_flds.taux, ocn.in_flds.tauy, ocn.τx, ocn.τy, direction=:Forward)

    if adv_scheme == :ekman_all_in_ML


        @loop_hor ocn i j let

            ϵ = ocn.ϵs[i, j]
            f = ocn.fs[i, j]

            τx = ocn.τx[i, j]
            τy = ocn.τy[i, j]

            h_ML = ocn.h_ML[i, j]
            Nz   = ocn.Nz[i, j] 
            s2ρh = ρ * h_ML * (ϵ^2.0 + f^2.0)

            ek_u = (ϵ * τx + f * τy) / s2ρh
            ek_v = (ϵ * τy - f * τx) / s2ρh

            FLDO = ocn.FLDO[i, j]

            if FLDO == -1
                ocn.u[:, i, j] .= ek_u
                ocn.v[:, i, j] .= ek_v
            else
                ocn.u[1:FLDO-1, i, j] .= ek_u
                ocn.v[1:FLDO-1, i, j] .= ek_v

                ocn.u[FLDO:Nz, i, j] .= 0.0
                ocn.v[FLDO:Nz, i, j] .= 0.0
            end
            
        end

    elseif adv_scheme == :ekman_simple_partition

        @loop_hor ocn i j let

            h_ML = ocn.h_ML[i, j]
#            s̃ = ocn.ϵs[i, j] + ocn.fs[i, j] * im
            s̃ = 1e-6 + ocn.fs[i, j] * im
            H̃ = √(1e-2 / s̃)
            H = abs(H̃)
            p̃ = exp(- ocn.h_ML[i, j] * H̃)
            
            M̃ = (ocn.τx[i, j] + ocn.τy[i, j] * im) / (ρ * s̃)
            M̃_DO = M̃ * p̃

            ṽ_ML = (M̃ - M̃_DO) / h_ML

            u_ML, v_ML = real(ṽ_ML), imag(ṽ_ML)

            FLDO = ocn.FLDO[i, j]



            if FLDO == -1
            
                ocn.u[:, i, j] .= u_ML
                ocn.v[:, i, j] .= v_ML

            else

                ocn.u[:, i, j] .= 0.0
                ocn.v[:, i, j] .= 0.0

                if FLDO > 1
                    ocn.u[1:FLDO-1, i, j] .= u_ML
                    ocn.v[1:FLDO-1, i, j] .= v_ML
                end

                H *= 3
                
                eklayer = getLayerFromDepth(
                    z  = - h_ML - H,
                    zs = ocn.cols.zs[i, j],  
                    Nz = ocn.Nz[i, j],
                )
            
                ṽ_DO = (M̃ - M̃_DO) / H
                u_DO, v_DO = real(ṽ_DO), imag(ṽ_DO)

                Δh     = ocn.hs[FLDO, i, j]
                Δh_top = h_ML + ocn.zs[FLDO, i, j]
                Δh_bot = Δh - Δh_top

                ocn.u[FLDO, i, j] = (Δh_top * u_ML + Δh_bot * u_DO) / Δh
                ocn.v[FLDO, i, j] = (Δh_top * v_ML + Δh_bot * v_DO) / Δh

                if FLDO != ocn.Nz[i, j]

                    if eklayer == -1
                        ocn.u[FLDO+1:end, i, j] .= u_DO
                        ocn.v[FLDO+1:end, i, j] .= v_DO
                    else
                        ocn.u[FLDO+1:eklayer, i, j] .= u_DO
                        ocn.v[FLDO+1:eklayer, i, j] .= v_DO
                    end

                end
                


            end
            
        end

    end
        
    # Calculate ∇⋅v
    for k=1:ocn.Nz_bone
        DisplacedPoleCoordinate.DIV!(ocn.gi, ocn.lays.u[k],  ocn.lays.v[k],  ocn.lays.div[k], ocn.lays.mask3[k])
    end

    # Calculate w
    for grid_idx in 1:size(ocn.valid_idx)[2]

        i = ocn.valid_idx[1, grid_idx]
        j = ocn.valid_idx[2, grid_idx]

        Nz = ocn.Nz[i, j]
        ocn.w[1, i, j] = 0.0
        for k = 2:Nz+1
            ocn.w[k, i, j] = ocn.w[k-1, i, j] + ocn.div[k-1, i, j]
        end
    end



    # Determine the temperature / salinity of FLDO layer
    for grid_idx in 1:size(ocn.valid_idx)[2]

        i = ocn.valid_idx[1, grid_idx]
        j = ocn.valid_idx[2, grid_idx]

        FLDO = ocn.FLDO[i, j]

        if FLDO != -1

            Δh     = ocn.hs[FLDO, i, j]
            Δh_top = ocn.h_ML[i, j] + ocn.zs[FLDO, i, j]
            Δh_bot = Δh - Δh_top

            ocn.Ts[FLDO, i, j] =  (Δh_top * ocn.T_ML[i, j] + Δh_bot * ocn.Ts[FLDO, i, j]) / Δh
            ocn.Ss[FLDO, i, j] =  (Δh_top * ocn.S_ML[i, j] + Δh_bot * ocn.Ss[FLDO, i, j]) / Δh

        end

    end

end

function stepOcean_Flow!(
    ocn  :: Ocean;
    cfgs...
)
    
    adv_scheme = cfgs[:adv_scheme]
    Δt         = cfgs[:Δt]

    if adv_scheme == :static
        return
    end



    # Pseudo code
    # 1. assign velocity field
    # 2. calculate temperature & salinity flux
    # 3. calculate temperature & salinity flux divergence
    # Gov eqn adv + diff: ∂T/∂t = - 1 / (ρ H1) ( ∇⋅(M1 T1) - (∇⋅M1) Tmid )
  
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

    for grid_idx in 1:size(ocn.valid_idx)[2]

        i = ocn.valid_idx[1, grid_idx]
        j = ocn.valid_idx[2, grid_idx]

        vadv_upwind!(
            ocn.cols.T_vadvs[i, j],
            ocn.cols.w[i, j],
            ocn.cols.Ts[i, j],
            ocn.cols.Δzs[i, j],
            ocn.Nz[i, j],
        )

        vadv_upwind!(
            ocn.cols.S_vadvs[i, j],
            ocn.cols.w[i, j],
            ocn.cols.Ss[i, j],
            ocn.cols.Δzs[i, j],
            ocn.Nz[i, j],
        )


    end

    for grid_idx in 1:size(ocn.valid_idx)[2]

        i = ocn.valid_idx[1, grid_idx]
        j = ocn.valid_idx[2, grid_idx]

        for k = 1:ocn.Nz[i, j]
            ocn.Ts[k, i, j] += Δt * ( ocn.T_vadvs[k, i, j] + ocn.T_hadvs[k, i, j] )
            ocn.Ss[k, i, j] += Δt * ( ocn.S_vadvs[k, i, j] + ocn.S_hadvs[k, i, j] )
        end


        zs   = ocn.cols.zs[i, j]
        hs   = ocn.cols.hs[i, j]
        h_ML = ocn.h_ML[i, j]
        FLDO = ocn.FLDO[i, j]
        Nz   = ocn.Nz[i, j]

        # Remix top layers
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

        OC_updateB!(ocn, i, j)
        OC_doConvectiveAdjustment!(ocn, i, j)


    end


end



function vadv_upwind!(
    vadvs  :: AbstractArray{Float64, 1},
    ws     :: AbstractArray{Float64, 1},
    qs     :: AbstractArray{Float64, 1},
    Δzs    :: AbstractArray{Float64, 1},
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
        return
    end

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

    # Do not update final layer
    vadvs[Nz] = 0.0 
    
end


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

#=
"""

Calculation follows Lax-Wendroff scheme.

Reference website:
  https://en.wikipedia.org/wiki/Lax%E2%80%93Wendroff_method

"""
function doZAdvection_MacCormack!(;
    Nz     :: Integer,
    qs     :: AbstractArray{Float64, 1},
    ws     :: AbstractArray{Float64, 1},
    Δhs    :: AbstractArray{Float64, 1},
    Δzs    :: AbstractArray{Float64, 1},
    Δt     :: Float64,
    qstmp  :: AbstractArray{Float64, 1},
    flxtmp :: AbstractArray{Float64, 1},
)

    #
    # Array information:
    # 
    # length(ws)     == length(qs) or Nz    + 1
    # length(Δhs)    == length(qs) or Nz
    # length(Δzs)    == length(qs) or Nz    - 1
    # length(qstmp)  == length(qs) or Nz    + 1
    # length(flxtmp) == length(qs) or Nz    + 1
    # 
    # Nz reveals if there is bottom of ocean
    # 

    # Step 1: calculate flux
    for k = 1:Nz
        flxtmp[k] = (ws[k] + ws[k+1]) / 2.0  * qs[k]
    end

    # Step 2: Lax step. Calculate staggered qs
    qstmp[1] = qs[1] - Δt / Δzs[
    for k = 2:length(qstmp)
        qstmp[k] = 
    end
   
    # Step 2: flux from bottom
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
    qs[1] = (qstmp[1] + qs[1]) / 2.0
    for k = 2:length(qstmp)-1
        qs[k] = (qstmp[k] + qs[k]) / 2.0 - Δt / Δzs[k-1] / 2.0 * (flxtmp[k-1] - flxtmp[k])
    end
   
    # Don't update qs[end] 
    
end
=#


