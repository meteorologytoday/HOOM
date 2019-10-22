


function calDiffAdv_QUICK!(
    ocn :: Ocean,
)

    calGRAD_CURV!(
        gi         = ocn.gi,
        Nx         = ocn.Nx,
        Ny         = ocn.Ny,
        Nz         = ocn.Nz,
        qs         = ocn.Ts,
        GRAD_bnd_x = ocn.GRAD_bnd_x,
        GRAD_bnd_y = ocn.GRAD_bnd_y,
        GRAD_bnd_z = ocn.GRAD_bnd_z,
        CURV_x     = ocn.CURV_x,
        CURV_y     = ocn.CURV_y,
        CURV_z     = ocn.CURV_z,
        mask3      = ocn.mask3,
        Δzs        = ocn.Δzs,
        hs         = ocn.hs,
    )

    calFluxDensity!(
        gi         = ocn.gi,
        Nx         = ocn.Nx,
        Ny         = ocn.Ny,
        Nz         = ocn.Nz,
        wq_bnd     = ocn.wT,
        qs         = ocn.Ts,
        GRAD_bnd_x = ocn.GRAD_bnd_x,
        GRAD_bnd_y = ocn.GRAD_bnd_y,
        GRAD_bnd_z = ocn.GRAD_bnd_z,
        CURV_x     = ocn.CURV_x,
        CURV_y     = ocn.CURV_y,
        CURV_z     = ocn.CURV_z,
        FLUX_DEN_x = ocn.FLUX_DEN_x,
        FLUX_DEN_y = ocn.FLUX_DEN_y,
        FLUX_DEN_z = ocn.FLUX_DEN_z,
        u_bnd      = ocn.u_bnd,
        v_bnd      = ocn.v_bnd,
        w_bnd      = ocn.w_bnd,
        mask3      = ocn.mask3,
        Δzs        = ocn.Δzs,
        D_hor      = ocn.Dh_T,
        D_ver      = ocn.Dv_T,
    )

    #=
    if any(isnan.(ocn.FLUX_DEN_x))
        throw(ErrorException("FLUX_DEN_x NaN"))
    end

    if any(isnan.(ocn.FLUX_DEN_y))
        throw(ErrorException("FLUX_DEN_y NaN"))
    end


    if any(isnan.(ocn.FLUX_DEN_z))
        throw(ErrorException("FLUX_DEN_z NaN"))
    end
    =#


    calTotalChange!(
        total_chg  = ocn.T_hadvs,
        gi         = ocn.gi,
        Nx         = ocn.Nx,
        Ny         = ocn.Ny,
        Nz         = ocn.Nz,
        FLUX_DEN_x = ocn.FLUX_DEN_x,
        FLUX_DEN_y = ocn.FLUX_DEN_y,
        FLUX_DEN_z = ocn.FLUX_DEN_z,
        mask3      = ocn.mask3,
        hs         = ocn.hs,
    )
end


function calTotalChange!(;
    total_chg  :: AbstractArray{Float64, 3},     # ( Nz_bone,  , Nx   , Ny   )
    gi         :: DisplacedPoleCoordinate.GridInfo,
    Nx         :: Integer,
    Ny         :: Integer,
    Nz         :: AbstractArray{Int64, 2},
    FLUX_DEN_x :: AbstractArray{Float64, 3},     # ( Nz_bone   ,  Nx+1, Ny   )
    FLUX_DEN_y :: AbstractArray{Float64, 3},     # ( Nz_bone   ,  Nx  , Ny+1 )
    FLUX_DEN_z :: AbstractArray{Float64, 3},     # ( Nz_bone+1 ,  Nx  , Ny   )
    mask3      :: AbstractArray{Float64, 3},     # ( Nz_bone   ,  Nx  , Ny   )
    hs         :: AbstractArray{Float64, 3},     # ( Nz_bone   ,  Nx  , Ny   )
)


    tmp = 0.0
    tmp_wT = 0.0
    tmp_v = 0.0
    tmp_σ = 0.0

    for i=1:Nx, j=1:Ny
        
        if mask3[1, i, j] == 0
            continue
        end

        i_e = (i==Nx) ? 1 : i+1


        if FLUX_DEN_z[1, i, j] != 0
            println("i: ", i, "; j:", j)
            throw(ErrorException("FLUX_DEN_z != 0.0"))
        end

        tmp_wT += FLUX_DEN_z[Nz[i, j]+1, i, j] * gi.dσ[i, j]
        tmp_σ += gi.dσ[i, j]

        for k=1:Nz[i, j]
            total_chg[k, i, j] = (
                - (
                     FLUX_DEN_x[k, i+1, j] * gi.ds4[i_e, j] - FLUX_DEN_x[k, i, j] * gi.ds4[i, j]
                   + FLUX_DEN_y[k, i, j+1] * gi.ds1[i, j+1] - FLUX_DEN_y[k, i, j] * gi.ds1[i, j]
                ) / gi.dσ[i, j]
                - (
                     FLUX_DEN_z[k, i, j] - FLUX_DEN_z[k+1, i, j]
                ) / hs[k, i, j]
            )

           if i==1 && FLUX_DEN_x[k, 1, j] != FLUX_DEN_x[k, Nx+1, j]
                println("i: ", i, "; j:", j)
                throw(ErrorException("FLUX_DEN_x does not match"))
           end
 
           if j==1 && ( FLUX_DEN_y[k, i, 1] != 0 ||  FLUX_DEN_x[k, i, Ny+1] != 0)
                println("i: ", i, "; j:", j)
                throw(ErrorException("FLUX_DEN_y != 0"))
           end
 

#=
            if i < Nx-1 && gi.ds2[i, j] != gi.ds4[i+1, j]
                println("i: ", i, "; j:", j)
                throw(ErrorException("ds2 ds4 does not match"))
            end
            
            if j < Ny-1 && gi.ds3[i, j] != gi.ds1[i, j+1]
                println("i: ", i, "; j:", j)
                println("ds3: ", gi.ds3[i, j], "; ds1: ", gi.ds1[i, j+1])
                throw(ErrorException("ds1 ds3 does not match"))
            end
=#

            tmp += total_chg[k, i, j] * hs[k, i, j] * gi.dσ[i, j]
            tmp_v += hs[k, i, j] * gi.dσ[i, j]

#            if (k, j) == (2, 10)
#                println(FLUX_DEN_x[k, i, j] * gi.ds4[i_e, j], " ::: ", FLUX_DEN_x[k, i+1, j] * gi.ds4[i, j])
#                tmp += FLUX_DEN_x[k, i+1, j] * gi.ds4[i_e, j] - FLUX_DEN_x[k, i, j] * gi.ds4[i, j]
#            end
            
#            if (k, i) == (1, 10) 
#                println(FLUX_DEN_y[k, i, j+1] * gi.ds1[i, j+1], " ::: ", FLUX_DEN_y[k, i, j] * gi.ds1[i, j])
#            end

        end

    end

#    println("SUM of total_chg weighted by volume: ", tmp, " / ", tmp_v, " = ", tmp/tmp_v)
#    println("wT total: ", tmp_wT)
    println("tmp_wT * ρc = ", tmp_wT * ρc / tmp_σ,"; If consider the affect of wT: ", (tmp - tmp_wT) /tmp_v)
end


function calFluxDensity!(;
    gi         :: DisplacedPoleCoordinate.GridInfo,
    Nx         :: Integer,
    Ny         :: Integer,
    Nz         :: AbstractArray{Int64, 2},
    wq_bnd     :: AbstractArray{Float64, 2},     # ( Nx  , Ny   )
    qs         :: AbstractArray{Float64, 3},     # ( Nz_bone   ,  Nx  , Ny   )
    u_bnd      :: AbstractArray{Float64, 3},     # ( Nz_bone   ,  Nx+1, Ny   )
    v_bnd      :: AbstractArray{Float64, 3},     # ( Nz_bone   ,  Nx  , Ny+1 )
    w_bnd      :: AbstractArray{Float64, 3},     # ( Nz_bone+1 ,  Nx  , Ny   )
    GRAD_bnd_x :: AbstractArray{Float64, 3},     # ( Nz_bone   ,  Nx+1, Ny   )
    GRAD_bnd_y :: AbstractArray{Float64, 3},     # ( Nz_bone   ,  Nx  , Ny+1 )
    GRAD_bnd_z :: AbstractArray{Float64, 3},     # ( Nz_bone+1 ,  Nx  , Ny   )
    CURV_x     :: AbstractArray{Float64, 3},     # ( Nz_bone   ,  Nx  , Ny   )
    CURV_y     :: AbstractArray{Float64, 3},     # ( Nz_bone   ,  Nx  , Ny   )
    CURV_z     :: AbstractArray{Float64, 3},     # ( Nz_bone   ,  Nx  , Ny   )
    FLUX_DEN_x :: AbstractArray{Float64, 3},     # ( Nz_bone   ,  Nx+1, Ny   )
    FLUX_DEN_y :: AbstractArray{Float64, 3},     # ( Nz_bone   ,  Nx  , Ny+1 )
    FLUX_DEN_z :: AbstractArray{Float64, 3},     # ( Nz_bone+1 ,  Nx  , Ny   )
    mask3      :: AbstractArray{Float64, 3},     # ( Nz_bone   ,  Nx  , Ny   )
    Δzs        :: AbstractArray{Float64, 3},     # ( Nz_bone-1 ,  Nx  , Ny   )
    D_hor      :: Float64,
    D_ver      :: Float64,
)

    # x
    for i=2:Nx, j=1:Ny 
        for k=1:Nz[i, j]
            if mask3[k, i, j] == 0.0 || mask3[k, i-1, j] == 0.0
                FLUX_DEN_x[k, i, j] = 0.0
            else
                q_star = (qs[k, i-1, j] + qs[k, i, j]) / 2.0 - gi.dx_w[i, j]^2.0/8.0 * ( ( u_bnd[k, i, j] >= 0.0 ) ? CURV_x[k, i-1, j] : CURV_x[k, i, j] )
                FLUX_DEN_x[k, i, j] = u_bnd[k, i, j] * q_star - D_hor * GRAD_bnd_x[k, i, j]
            end
        end
    end

    # x - periodic
    for j=1:Ny 
        for k=1:Nz[1, j]
            if mask3[k, 1, j] == 0.0 || mask3[k, Nx, j] == 0.0
                FLUX_DEN_x[k, 1, j] = FLUX_DEN_x[k, Nx+1, j] = 0.0
            else
                q_star = (qs[k, Nx, j] + qs[k, 1, j]) / 2.0 - gi.dx_w[1, j]^2.0/8.0 * ( ( u_bnd[k, 1, j] >= 0.0 ) ? CURV_x[k, Nx, j] : CURV_x[k, 1, j] )
                FLUX_DEN_x[k, 1, j] = FLUX_DEN_x[k, Nx+1, j] = u_bnd[k, 1, j] * q_star - D_hor * GRAD_bnd_x[k, 1, j]
            end
        end
    end


    # y
    for i=1:Nx, j=2:Ny-1
        for k=1:Nz[i, j]
            if mask3[k, i, j] == 0.0 || mask3[k, i, j-1] == 0.0
                FLUX_DEN_y[k, i, j] = 0.0
            else
                q_star = (qs[k, i, j-1] + qs[k, i, j]) / 2.0 - gi.dy_s[i, j]^2.0/8.0 * ( ( v_bnd[k, i, j] >= 0.0 ) ? CURV_y[k, i, j-1] : CURV_y[k, i, j] )
                FLUX_DEN_y[k, i, j] = v_bnd[k, i, j] * q_star - D_hor * GRAD_bnd_y[k, i, j]
            end
        end
    end


    # z
    for i=1:Nx, j=1:Ny

        if mask3[1, i, j] == 0.0
            continue
        end

        FLUX_DEN_z[1, i, j] = 0.0

        _Nz = Nz[i, j]

        for k=2:_Nz
            q_star = (qs[k, i, j] + qs[k-1, i, j]) / 2.0 - Δzs[k-1, i, j]^2.0/8.0 * ( ( w_bnd[k, i, j] >= 0.0 ) ? CURV_z[k, i, j] : CURV_z[k-1, i, j] )
            FLUX_DEN_z[k, i, j] = w_bnd[k, i, j] * q_star - D_ver * GRAD_bnd_z[k, i, j]
 
        end

        FLUX_DEN_z[_Nz+1, i, j] = wq_bnd[i, j] = FLUX_DEN_z[_Nz, i, j]
    end

end



function calGRAD_CURV!(;
    gi         :: DisplacedPoleCoordinate.GridInfo,
    Nx         :: Integer,
    Ny         :: Integer,
    Nz         :: AbstractArray{Int64, 2},
    qs         :: AbstractArray{Float64, 3},     # ( Nz_bone   ,  Nx  , Ny   )
    GRAD_bnd_x :: AbstractArray{Float64, 3},     # ( Nz_bone   ,  Nx+1, Ny   )
    GRAD_bnd_y :: AbstractArray{Float64, 3},     # ( Nz_bone   ,  Nx  , Ny+1 )
    GRAD_bnd_z :: AbstractArray{Float64, 3},     # ( Nz_bone+1 ,  Nx  , Ny   )
    CURV_x     :: AbstractArray{Float64, 3},     # ( Nz_bone   ,  Nx  , Ny   )
    CURV_y     :: AbstractArray{Float64, 3},     # ( Nz_bone   ,  Nx  , Ny   )
    CURV_z     :: AbstractArray{Float64, 3},     # ( Nz_bone   ,  Nx  , Ny   )
    mask3      :: AbstractArray{Float64, 3},     # ( Nz_bone   ,  Nx  , Ny   )
    Δzs        :: AbstractArray{Float64, 3},     # ( Nz_bone-1 ,  Nx  , Ny   )
    hs         :: AbstractArray{Float64, 3},     # ( Nz_bone   ,  Nx  , Ny   )
)

    # x 
    for i=2:Nx, j=1:Ny 
        for k=1:Nz[i, j]
            GRAD_bnd_x[k, i, j] = (
                ( mask3[k, i, j] == 0.0 || mask3[k, i-1, j] == 0.0 )  
                ? 0.0 : ( qs[k, i, j] - qs[k, i-1, j] ) / gi.dx_w[i, j] 
            )
        end
    end

    # x - periodic
    for j=1:Ny
        for k=1:Nz[1, j]
            GRAD_bnd_x[k, 1, j] = GRAD_bnd_x[k, Nx+1, j] = (
                ( mask3[k, 1, j] == 0.0 || mask3[k, Nx, j] == 0.0 )
                ? 0.0 : ( qs[k, 1, j] - qs[k, Nx, j] ) / gi.dx_w[1, j]
            )
        end
    end

    # y
    for i=1:Nx, j=2:Ny-1
        for k=1:Nz[i, j]
            GRAD_bnd_y[k, i, j] = (
                ( mask3[k, i, j] == 0.0 || mask3[k, i, j-1] == 0.0 )
                ? 0.0 : ( qs[k, i, j] - qs[k, i, j-1] ) / gi.dy_s[i, j]
            )
        end
    end

    # z
    for i=1:Nx, j=1:Ny

        if mask3[1, i, j] == 0.0
            continue
        end

        _Nz = Nz[i, j]
        GRAD_bnd_z[1, i, j] = GRAD_bnd_z[_Nz+1, i, j] = 0.0
        for k=2:_Nz
            GRAD_bnd_z[k, i, j] = ( qs[k-1, i, j] - qs[k, i, j] ) / Δzs[k-1, i, j]
        end

    end

    # CURV
    for i=1:Nx, j=1:Ny
        for k=1:Nz[i, j]
            CURV_x[k, i, j] = ( GRAD_bnd_x[k, i+1, j  ] - GRAD_bnd_x[k  , i, j] ) / gi.dx_c[i, j]
            CURV_y[k, i, j] = ( GRAD_bnd_y[k, i  , j+1] - GRAD_bnd_y[k  , i, j] ) / gi.dy_c[i, j]
            CURV_z[k, i, j] = ( GRAD_bnd_z[k, i  , j  ] - GRAD_bnd_z[k+1, i, j] ) / hs[k, i, j]
        end
    end

    #=
    if any(isnan.(GRAD_bnd_x))
        throw(ErrorException("GRAD_bnd_x NaN"))
    end

    if any(isnan.(GRAD_bnd_y))
        throw(ErrorException("GRAD_bnd_y NaN"))
    end

    if any(isnan.(GRAD_bnd_z))
        throw(ErrorException("GRAD_bnd_z NaN"))
    end

    =#

end



function calVerVelBnd!(;
    gi       :: DisplacedPoleCoordinate.GridInfo,
    Nx       :: Integer,
    Ny       :: Integer,
    Nz       :: AbstractArray{Int64, 2},
    u_bnd    :: AbstractArray{Float64, 3},   # ( Nz_bone  , Nx+1, Ny   )
    v_bnd    :: AbstractArray{Float64, 3},   # ( Nz_bone  , Nx  , Ny+1 )
    w_bnd    :: AbstractArray{Float64, 3},   # ( Nz_bone+1, Nx  , Ny   )
    hs       :: AbstractArray{Float64, 3},   # ( Nz_bone  , Nx  , Ny   )
    div      :: AbstractArray{Float64, 3},   # ( Nz_bone  , Nx  , Ny   )
    mask3    :: AbstractArray{Float64, 3},   # ( Nz_bone  , Nx  , Ny   )
)

    for i=1:Nx, j=1:Ny

        w_bnd[1, i, j] = 0.0

        for k=1:Nz[i, j]

            if mask3[k, i, j] == 0.0
                break
            end
            
            div[k, i, j] =  (  
                u_bnd[k, i+1, j  ]  * gi.ds2[i, j]
              - u_bnd[k, i,   j  ]  * gi.ds4[i, j]
              + v_bnd[k, i,   j+1]  * gi.ds3[i, j]
              - v_bnd[k, i,   j  ]  * gi.ds1[i, j]
            ) / gi.dσ[i, j]

            w_bnd[k+1, i, j] = w_bnd[k, i, j] + div[k, i, j] * hs[k, i, j]
        end
    end
 
end



function calHorVelBnd!(;
    Nx       :: Integer,
    Ny       :: Integer,
    Nz       :: AbstractArray{Int64, 2},
    weight_e :: AbstractArray{Float64, 2},   # (Nx+1, Ny)
    weight_n :: AbstractArray{Float64, 2},   # (Nx, Ny+1)
    u        :: AbstractArray{Float64, 3},   # (Nz_bone, Nx, Ny)
    v        :: AbstractArray{Float64, 3},   # (Nz_bone, Nx, Ny)
    u_bnd    :: AbstractArray{Float64, 3},   # (Nz_bone, Nx+1, Ny)
    v_bnd    :: AbstractArray{Float64, 3},   # (Nz_bone, Nx, Ny+1)
    mask3    :: AbstractArray{Float64, 3},   # (Nz_bone, Nx, Ny)
)

    # x
    for i=2:Nx, j=1:Ny
        for k=1:Nz[i, j]
            if mask3[k, i, j] == 0.0 || mask3[k, i-1, j] == 0.0
                u_bnd[k, i, j] = 0.0
            else
                u_bnd[k, i, j] = u[k, i-1, j] * (1.0 - weight_e[i, j]) + u[k, i, j] * weight_e[i, j]
            end
        end
    end
    
    # x - periodic
    for j=1:Ny
        for k=1:Nz[1, j]
            if mask3[k, 1, j] == 0.0 || mask3[k, Nx, j] == 0.0
                u_bnd[k, 1, j] = u_bnd[k, Nx+1, j] = 0.0
            else
                u_bnd[k, 1, j] = u_bnd[k, Nx+1, j] = u[k, Nx, j] * (1.0 - weight_e[1, j]) + u[k, 1, j] * weight_e[1, j]
            end
        end
    end

    # y
    for i=1:Nx, j=2:Ny-1
        for k=1:Nz[i, j]
            if mask3[k, i, j-1] == 0.0 || mask3[k, i, j] == 0.0
                v_bnd[k, i, j] = 0.0
            else
                v_bnd[k, i, j] = v[k, i, j-1] * (1.0 - weight_n[i, j]) + v[k, i, j] * weight_n[i, j]
            end
        end
    end

end
