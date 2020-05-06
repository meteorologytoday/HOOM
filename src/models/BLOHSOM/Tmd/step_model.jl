function determineVelocity!(
    m :: TmdModel,
)
    @fast_extract m

    wksp  = co.wksp

    div   = getSpace!(wksp, :T)
    tmp_T = getSpace!(wksp, :T)

    u_U = getSpace!(wksp, :U)
    v_V = getSpace!(wksp, :V)

    #=
    fr.u_U .= 0.0
    fr.v_V .= 0.0
    fr.u_U[1:10, :, 41:50] .= 2.0
=#
    # filter velocity
    mul_autoflat!(u_U, co.ASUM.filter_U, fr.u_U)
    mul_autoflat!(v_V, co.ASUM.filter_V, fr.v_V)

    fr.u_U .= u_U
    fr.v_V .= v_V
 
    calDIV!(
        ASUM = co.ASUM,
        u_U = fr.u_U,
        v_V = fr.v_V,
        div   = div,
        workspace = tmp_T
    )

    calVerVelBnd!(
        gi    = ev.gi,
        Nx    = ev.Nx,
        Ny    = ev.Ny,
        Nz    = ev.Nz_av,
        w_bnd = fr.w_W,
        hs    = co.Δz_T,
        div   = div,
        mask3 = ev.mask3,
    )
    

end


function advectTracer!(
    m   :: TmdModel,
)

    @fast_extract m

#    println("BEFORE:T[1:2]", st.X[1:2, 3, 3, 1], "; FDLO = ", st.FLDO[3,3], "; FLDO_ratio_top")


    calFLDOPartition!(m)

    for x = 1:ev.NX

        ΔX   = view(st.ΔX,      :, :, x)
        X_ML = view(st.X_ML,    :, :, x)
        X    = view(co.cols.X,  :, :, x)

        for I in CartesianIndices(X)
            ΔX[I] = mixFLDO!(
                qs   = X[I],
                zs   = co.cols.z_bnd_av[I],
                hs   = co.cols.Δz_T[I],
                q_ML = X_ML[I],
                FLDO = st.FLDO[I],
                FLDO_ratio_top = st.FLDO_ratio_top[I],
                FLDO_ratio_bot = st.FLDO_ratio_bot[I],
            )
        end
    end


#    println("AFTER:T[1:2]", st.X[1:2, 3, 3, 1], ";h_ML=", st.h_ML[3, 3])

    Δt = ev.Δt_substep
    tmp_T = getSpace!(co.wksp, :T)
    
    # Pseudo code
    # 1. calculate tracer flux
    # 2. calculate tracer flux divergence

    T_sum_old = sum(co.ASUM.T_Δvol_T * view(st.T, :))
    calDiffAdv_QUICKEST_SpeedUp!(m, Δt)


    @. st.X  += Δt * co.XFLUX_CONV
    @. st.ΔX += Δt * st.dΔXdt

#    T_sum_new = sum(co.ASUM.T_Δvol_T * view(st.T, :))

#    println("Tsum: ", T_sum_old, " => ", T_sum_new, "; change = ", (T_sum_new - T_sum_old) / T_sum_old * 100.0, " %")


#    println("## Before T_ML=", st.X_ML[3,3,1])
    for x = 1:ev.NX

        ΔX   = view(st.ΔX,      :, :, x)
        X_ML = view(st.X_ML,    :, :, x)
        X    = view(co.cols.X,  :, :, x)

        for I in CartesianIndices(X)
            X_ML[I] = unmixFLDOKeepDiff!(
                qs   = X[I],
                zs   = co.cols.z_bnd_av[I],
                hs   = co.cols.Δz_T[I],
                h_ML = st.h_ML[I],
                FLDO = st.FLDO[I],
                Nz   = ev.Nz_av[I],
                Δq   = ΔX[I],
            )
        end
    end

#    println("## After T_ML=", st.X_ML[3,3,1])

end

function calVerVelBnd!(;
    gi       :: PolelikeCoordinate.GridInfo,
    Nx       :: Integer,
    Ny       :: Integer,
    Nz       :: AbstractArray{Int64, 2},
    w_bnd    :: AbstractArray{Float64, 3},   # ( Nz_bone+1, Nx  , Ny   )
    hs       :: AbstractArray{Float64, 3},   # ( Nz_bone  , Nx  , Ny   )
    div      :: AbstractArray{Float64, 3},   # ( Nz_bone  , Nx  , Ny   )
    mask3    :: AbstractArray{Float64, 3},   # ( Nz_bone  , Nx  , Ny   )
)

#    local tmp = tmp_σ = 0.0
    for i=1:Nx, j=1:Ny
        
        _Nz = Nz[i, j]
        w_bnd[1, i, j]     = 0.0
        w_bnd[_Nz+1, i, j] = 0.0

        for k=_Nz:-1:2

            if mask3[k, i, j] == 0.0
                break
            end
            
            #w_bnd[k+1, i, j] = w_bnd[k, i, j] + div[k, i, j] * hs[k, i, j]
            w_bnd[k, i, j] = w_bnd[k+1, i, j] - div[k, i, j] * hs[k, i, j]
        end

#        tmp   += w_bnd[Nz[i, j]+1, i, j] * gi.dσ[i, j]
#        tmp_σ += gi.dσ[i, j]
    end

#   #println("tmp: ", tmp, "; tmp_σ: ", tmp_σ, "; Average w: ", tmp/tmp_σ)

end


function calMixedLayer_dΔqdt!(;
    Nx          :: Integer,
    Ny          :: Integer,
    Nz          :: AbstractArray{Int64, 2},
    FLUX_CONV_h :: AbstractArray{Float64, 3},     # ( Nz_bone  ,  Nx, Ny )
    FLUX_DEN_z  :: AbstractArray{Float64, 3},     # ( Nz_bone+1,  Nx, Ny )
    dΔqdt       :: AbstractArray{Float64, 2},     # ( Nx, Ny )
    mask        :: AbstractArray{Float64, 2},     # ( Nx, Ny )
    FLDO        :: AbstractArray{Int64, 2},       # ( Nx, Ny )
    h_ML        :: AbstractArray{Float64, 2},     # ( Nx, Ny )
    hs          :: AbstractArray{Float64, 3},     # ( Nz_bone  ,  Nx, Ny )
    zs          :: AbstractArray{Float64, 3},     # ( Nz_bone+1,  Nx, Ny )
) 

    for i=1:Nx, j=1:Ny

        if mask[i, j] == 0.0
            continue
        end

        _FLDO = FLDO[i, j]

        if _FLDO == -1
            continue
        end

        tmp = 0.0
        for k = 1:_FLDO-1
            tmp += FLUX_CONV_h[k, i, j] * hs[k, i, j]
        end
        tmp += ( 
              FLUX_CONV_h[_FLDO, i, j] * zs[_FLDO, i, j] 
            + ( FLUX_DEN_z[_FLDO+1, i, j] * zs[_FLDO, i, j] - FLUX_DEN_z[_FLDO, i, j] * zs[_FLDO+1, i, j] ) / hs[_FLDO, i, j]
        )

        dΔqdt[i, j] = tmp / h_ML[i, j]
    end

end


