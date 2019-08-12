using Statistics

function stepOcean_hz!(
    ocn  :: Ocean;
    cfgs...
)
 
    # Pseudo code
    # 1. assign velocity field
    # 2. calculate temperature & salinity flux
    # 3. calculate temperature & salinity flux divergence
    # Gov eqn adv + diff: ∂T/∂t = - 1 / (ρ H1) ( ∇⋅(M1 T1) - (∇⋅M1) Tmid )
   
    # Transform input wind stress vector first
    DisplacedPoleCoordinate.project!(ocn.gi, ocn.in_flds.taux, ocn.in_flds.tauy, ocn.τx, ocn.τy, direction=:Forward)



    for grid_idx in 1:size(ocn.valid_idx)[2]

        i = ocn.valid_idx[1, grid_idx]
        j = ocn.valid_idx[2, grid_idx]

        ϵ = 1e-6 #ocn.ϵs[i, j]
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
        
        #for k=1:Nz
        #    ocn.uT[k, i, j] = ocn.u[k, i, j] * ocn.Ts[k, i, j]
        #    ocn.vT[k, i, j] = ocn.v[k, i, j] * ocn.Ts[k, i, j]
        #    ocn.uS[k, i, j] = ocn.u[k, i, j] * ocn.Ss[k, i, j]
        #    ocn.vS[k, i, j] = ocn.v[k, i, j] * ocn.Ss[k, i, j]
        #end


    end
        
    #println(sum(ocn.u[isfinite.(ocn.u)])) 



    for k=1:ocn.Nz_bone
        # Calculate ∇⋅v
        DisplacedPoleCoordinate.DIV!(ocn.gi, ocn.lays.u[k],  ocn.lays.v[k],  ocn.lays.div[k])
        
        # Calculate ∇⋅(vT)
        DisplacedPoleCoordinate.DIV!(ocn.gi, ocn.lays.uT[k], ocn.lays.vT[k], ocn.lays.divTflx[k])

        # Calculate ∇⋅(vS)
        DisplacedPoleCoordinate.DIV!(ocn.gi, ocn.lays.uS[k], ocn.lays.vS[k], ocn.lays.divSflx[k])
        
        # Calculate ∇∇T, ∇∇S
        DisplacedPoleCoordinate.∇∇!(ocn.gi, ocn.lays.Ts[k], ocn.lays.∇∇T[k])
        DisplacedPoleCoordinate.∇∇!(ocn.gi, ocn.lays.Ss[k], ocn.lays.∇∇S[k])
    end


    for grid_idx in 1:size(ocn.valid_idx)[2]

        i = ocn.valid_idx[1, grid_idx]
        j = ocn.valid_idx[2, grid_idx]

        Nz = ocn.Nz[i, j]
        ocn.w[1, i, j] = 0.0
        for k = 2:Nz+1
            ocn.w[k, i, j] = ocn.w[k-1, i, j] + ocn.div[k-1, i, j]
        end
    end


    # Now doing advection
    
end
