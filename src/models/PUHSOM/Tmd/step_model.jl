function advectTracer!(
    model   :: TmdModel,
)


    Δt = model.env.Δt

    state = model.state
    core  = model.core
    env   = model.env
    wksp  = core.wksp

    div   = getSpace!(wksp, :T)
    tmp_T = getSpace!(wksp, :T)

    calDIV!(
        ASUM = core.ASUM,
        u_bnd = state.u_U,
        v_bnd = state.v_V,
        div   = div,
        workspace = tmp_T
    )

    calVerVelBnd!(
        gi    = env.gi,
        Nx    = env.Nx,
        Ny    = env.Ny,
        Nz    = env.Nz_av,
        w_bnd = state.w_W,
        hs    = core.dz_W,
        div   = div,
        mask3 = env.mask3,
    )
   

    # Pseudo code
    # 1. calculate tracer flux
    # 2. calculate tracer flux divergence
    calDiffAdv_QUICKEST_SpeedUp!(model, Δt)

#=
    for x=1:env.NX, i = 1:env.Nx, j=1:env.Ny
            for k = 1:env.Nz_av[i, j]
                state.X[k, i, j, x] += Δt * core.XFLUX_CONV[k, i, j, x]
            end
    end
=#
    let
        for x=1:env.NX
            tmp_T .= view(core.XFLUX_CONV, :, :, :, x)
            tmp_T .*= Δt
            state.X[:, :, :, x] .+= tmp_T
        end
    end

end

function calDiffAdv_QUICKEST_SpeedUp!(
    model       :: TmdModel,
    Δt          :: Float64,
)
   
    core = model.core
    env     = model.env
    state   = model.state
    ASUM    = core.ASUM
    wksp    = core.wksp

    GRAD_bnd_x = getSpace!(wksp, :U)
    GRAD_bnd_y = getSpace!(wksp, :V)
    GRAD_bnd_z = getSpace!(wksp, :W)
        
    CURV_x = getSpace!(wksp, :T)
    CURV_y = getSpace!(wksp, :T)
    CURV_z = getSpace!(wksp, :T)
        
    tmp1 = getSpace!(wksp, :T)
    tmp2 = getSpace!(wksp, :T)
    tmp3 = getSpace!(wksp, :T)

    for x=1:env.NX
 
        X            = view(model.state.X,     :, :, :, x)

        XFLUX_bot    = view(core.XFLUX_bot,       :, :, x)
        XFLUX_CONV   = view(core.XFLUX_CONV,   :, :, :, x)
        XFLUX_CONV_h = view(core.XFLUX_CONV_h, :, :, :, x)
        XFLUX_DEN_x  = view(core.XFLUX_DEN_x,  :, :, :, x)
        XFLUX_DEN_y  = view(core.XFLUX_DEN_y,  :, :, :, x)
        XFLUX_DEN_z  = view(core.XFLUX_DEN_z,  :, :, :, x)

        let
            mul!(view(GRAD_bnd_x, :), ASUM.mtx_GRAD_X, view(X, :))
            mul!(view(GRAD_bnd_y, :), ASUM.mtx_GRAD_Y, view(X, :))
            mul!(view(GRAD_bnd_z, :), ASUM.mtx_GRAD_Z, view(X, :))
     
            mul!(view(CURV_x, :), ASUM.mtx_CURV_X, view(GRAD_bnd_x, :))
            mul!(view(CURV_y, :), ASUM.mtx_CURV_Y, view(GRAD_bnd_y, :))
            mul!(view(CURV_z, :), ASUM.mtx_CURV_Z, view(GRAD_bnd_z, :))
     
        end

        #println("Flux Density")
        calFluxDensity!(
            gi         = env.gi,
            Nx         = env.Nx,
            Ny         = env.Ny,
            Nz         = env.Nz_av,
            FLUX_bot   = XFLUX_bot,
            qs         = X,
            GRAD_bnd_x = GRAD_bnd_x,
            GRAD_bnd_y = GRAD_bnd_y,
            GRAD_bnd_z = GRAD_bnd_z,
            CURV_x     = CURV_x,
            CURV_y     = CURV_y,
            CURV_z     = CURV_z,
            FLUX_DEN_x = XFLUX_DEN_x,
            FLUX_DEN_y = XFLUX_DEN_y,
            FLUX_DEN_z = XFLUX_DEN_z,
            u_bnd      = state.u_U,
            v_bnd      = state.v_V,
            w_bnd      = state.w_W,
            mask3          = env.mask3,
            noflux_x_mask3 = env.noflux_x_mask3,
            noflux_y_mask3 = env.noflux_y_mask3,
            Δzs        = core.dz_T,
            D_hor      = env.Kh_X[x],
            D_ver      = env.Kv_X[x],
            Δt         = Δt,
        )

    


   #println("TOTAL CHANGE")
        let
            mul!(view(XFLUX_CONV_h, :), ASUM.mtx_DIV_X, view(XFLUX_DEN_x, :))
            mul!(view(tmp1, :),   ASUM.mtx_DIV_Y, view(XFLUX_DEN_y, :))
            mul!(view(tmp2, :),   ASUM.mtx_DIV_Z, view(XFLUX_DEN_z, :))

            XFLUX_CONV_h .+= tmp1
            XFLUX_CONV_h .*= -1.0
            XFLUX_CONV .= XFLUX_CONV_h 
            XFLUX_CONV .-= tmp2

    #        for j=1:ocn.Ny, i=1:ocn.Nx, k=1:ocn.Nz_bone
    #            FLUX_CONV_h[k, i, j] = - ( ocn.workspace1[k, i, j] + ocn.workspace2[k, i, j] )
    #            FLUX_CONV[k, i, j] = FLUX_CONV_h[k, i, j] - ocn.workspace3[k, i, j]
    #        end
        end

    end

end


function calTotalChange!(;
    FLUX_CONV  :: AbstractArray{Float64, 3},     # ( Nz_bone,  , Nx   , Ny   )
    FLUX_CONV_h:: AbstractArray{Float64, 3},     # ( Nz_bone,  , Nx   , Ny   )
    gi         :: PolelikeCoordinate.GridInfo,
    Nx         :: Integer,
    Ny         :: Integer,
    Nz         :: AbstractArray{Int64, 2},
    FLUX_DEN_x :: AbstractArray{Float64, 3},     # ( Nz_bone   ,  Nx, Ny   )
    FLUX_DEN_y :: AbstractArray{Float64, 3},     # ( Nz_bone   ,  Nx  , Ny+1 )
    FLUX_DEN_z :: AbstractArray{Float64, 3},     # ( Nz_bone+1 ,  Nx  , Ny   )
    mask3      :: AbstractArray{Float64, 3},     # ( Nz_bone   ,  Nx  , Ny   )
    hs         :: AbstractArray{Float64, 3},     # ( Nz_bone   ,  Nx  , Ny   )
)


    #tmp = 0.0
    #tmp_wT = 0.0
    #tmp_v = 0.0
    #tmp_σ = 0.0

    for i=1:Nx, j=1:Ny
        
        if mask3[1, i, j] == 0
            continue
        end

        #i_e = (i==Nx) ? 1 : i+1


 #       if FLUX_DEN_z[1, i, j] != 0
 #          #println("i: ", i, "; j:", j)
 #           throw(ErrorException("FLUX_DEN_z != 0.0"))
 #       end

      # tmp_wT += FLUX_DEN_z[Nz[i, j]+1, i, j] * gi.dσ[i, j]
      # tmp_σ += gi.dσ[i, j]

#=
        if (i, j) == (48, 89)
           #println("FLUX_X = ", FLUX_DEN_x[1:6, i+1, j], "; ", FLUX_DEN_x[1:6, i, j])
           #println("FLUX_X conv=", FLUX_DEN_x[1:6, i+1, j] * gi.DY[i+1, j] - FLUX_DEN_x[1:6, i, j] * gi.DY[i, j])
           #println("FLUX_y conv=", FLUX_DEN_y[1:6, i, j+1] * gi.DX[i, j+1] - FLUX_DEN_y[1:6, i, j] * gi.DX[i, j])

           #println("DX: ", gi.DX[i, j:j+1])
           #println("DY: ", gi.DY[i:i+1, j])
        end
=#

        for k=1:Nz[i, j]

            _CONV_h = (
                - (
                     FLUX_DEN_x[k, i+1, j] * gi.DY[i+1, j] - FLUX_DEN_x[k, i, j] * gi.DY[i, j]
                   + FLUX_DEN_y[k, i, j+1] * gi.DX[i, j+1] - FLUX_DEN_y[k, i, j] * gi.DX[i, j]
                ) / gi.dσ[i, j]
            )

            FLUX_CONV_h[k, i, j] = _CONV_h

            FLUX_CONV[k, i, j] = ( 
                _CONV_h - (
                     FLUX_DEN_z[k, i, j] - FLUX_DEN_z[k+1, i, j]
                ) / hs[k, i, j]
            )

#=
           if i==1 && FLUX_DEN_x[k, 1, j] != FLUX_DEN_x[k, Nx+1, j]
               #println("i: ", i, "; j:", j)
                throw(ErrorException("FLUX_DEN_x does not match"))
           end
 
           if j==1 && ( FLUX_DEN_y[k, i, 1] != 0 ||  FLUX_DEN_x[k, i, Ny+1] != 0)
               #println("i: ", i, "; j:", j)
                throw(ErrorException("FLUX_DEN_y != 0"))
           end
 
=#
#=
            if i < Nx-1 && gi.ds2[i, j] != gi.ds4[i+1, j]
               #println("i: ", i, "; j:", j)
                throw(ErrorException("ds2 ds4 does not match"))
            end
            
            if j < Ny-1 && gi.ds3[i, j] != gi.ds1[i, j+1]
               #println("i: ", i, "; j:", j)
               #println("ds3: ", gi.ds3[i, j], "; ds1: ", gi.ds1[i, j+1])
                throw(ErrorException("ds1 ds3 does not match"))
            end
=#

     #       tmp += FLUX_CONV[k, i, j] * hs[k, i, j] * gi.dσ[i, j]
     #       tmp_v += hs[k, i, j] * gi.dσ[i, j]

#            if (k, j) == (2, 10)
#               #println(FLUX_DEN_x[k, i, j] * gi.ds4[i_e, j], " ::: ", FLUX_DEN_x[k, i+1, j] * gi.ds4[i, j])
#                tmp += FLUX_DEN_x[k, i+1, j] * gi.ds4[i_e, j] - FLUX_DEN_x[k, i, j] * gi.ds4[i, j]
#            end
            
#            if (k, i) == (1, 10) 
#               #println(FLUX_DEN_y[k, i, j+1] * gi.ds1[i, j+1], " ::: ", FLUX_DEN_y[k, i, j] * gi.ds1[i, j])
#            end

        end

    end

    #println("SUM of FLUX_CONV weighted by volume: ", tmp, " / ", tmp_v, " = ", tmp/tmp_v)
    #println("wQ total: ", tmp_wT/tmp_σ)
    #println("If consider the affect of wQ: ", (tmp - tmp_wT) /tmp_v)
end

function calFluxDensity!(;
    gi         :: PolelikeCoordinate.GridInfo,
    Nx         :: Integer,
    Ny         :: Integer,
    Nz         :: AbstractArray{Int64, 2},
    FLUX_bot     :: AbstractArray{Float64, 2},     # ( Nx  , Ny   )
    qs         :: AbstractArray{Float64, 3},     # ( Nz_bone   ,  Nx  , Ny   )
    u_bnd      :: AbstractArray{Float64, 3},     # ( Nz_bone   ,  Nx, Ny   )
    v_bnd      :: AbstractArray{Float64, 3},     # ( Nz_bone   ,  Nx  , Ny+1 )
    w_bnd      :: AbstractArray{Float64, 3},     # ( Nz_bone+1 ,  Nx  , Ny   )
    GRAD_bnd_x :: AbstractArray{Float64, 3},     # ( Nz_bone   ,  Nx, Ny   )
    GRAD_bnd_y :: AbstractArray{Float64, 3},     # ( Nz_bone   ,  Nx  , Ny+1 )
    GRAD_bnd_z :: AbstractArray{Float64, 3},     # ( Nz_bone+1 ,  Nx  , Ny   )
    CURV_x     :: AbstractArray{Float64, 3},     # ( Nz_bone   ,  Nx  , Ny   )
    CURV_y     :: AbstractArray{Float64, 3},     # ( Nz_bone   ,  Nx  , Ny   )
    CURV_z     :: AbstractArray{Float64, 3},     # ( Nz_bone   ,  Nx  , Ny   )
    FLUX_DEN_x :: AbstractArray{Float64, 3},     # ( Nz_bone   ,  Nx, Ny   )
    FLUX_DEN_y :: AbstractArray{Float64, 3},     # ( Nz_bone   ,  Nx  , Ny+1 )
    FLUX_DEN_z :: AbstractArray{Float64, 3},     # ( Nz_bone+1 ,  Nx  , Ny   )
    mask3      :: AbstractArray{Float64, 3},     # ( Nz_bone   ,  Nx  , Ny   )
    noflux_x_mask3 :: AbstractArray{Float64, 3}, # ( Nz_bone   ,  Nx, Ny   )
    noflux_y_mask3 :: AbstractArray{Float64, 3}, # ( Nz_bone   ,  Nx  , Ny+1 )
    Δzs        :: AbstractArray{Float64, 3},     # ( Nz_bone-1 ,  Nx  , Ny   )
    D_hor      :: Float64,
    D_ver      :: Float64,
    Δt         :: Float64,
)

    # x
    for i=2:Nx, j=1:Ny
        for k=1:Nz[i, j]
            if noflux_x_mask3[k, i, j] == 0.0
                FLUX_DEN_x[k, i, j] = 0.0
            else
                CURV_r = ( u_bnd[k, i, j] >= 0.0 ) ? CURV_x[k, i-1, j] : CURV_x[k, i, j]
                uΔt    = u_bnd[k, i, j] * Δt
                q_star = (qs[k, i-1, j] + qs[k, i, j]) / 2.0 - uΔt / 2.0 * GRAD_bnd_x[k, i, j] + ( D_hor * Δt / 2.0 - gi.dx_w[i, j]^2.0/6.0 + uΔt^2.0 / 6.0 ) * CURV_r

                FLUX_DEN_x[k, i, j] = u_bnd[k, i, j] * q_star - D_hor * ( GRAD_bnd_x[k, i, j] - uΔt / 2.0 * CURV_r )
            end
        end
    end

    # x - periodic
    for j=1:Ny 
        for k=1:Nz[1, j]
            if noflux_x_mask3[k, 1, j] == 0.0
                FLUX_DEN_x[k, 1, j] = 0.0
            else

                CURV_r = ( u_bnd[k, 1, j] >= 0.0 ) ? CURV_x[k, Nx, j] : CURV_x[k, 1, j]
                uΔt    = u_bnd[k, 1, j] * Δt
                q_star = (qs[k, Nx, j] + qs[k, 1, j]) / 2.0 - uΔt / 2.0 * GRAD_bnd_x[k, 1, j] + ( D_hor * Δt / 2.0 - gi.dx_w[1, j]^2.0/6.0 + uΔt^2.0 / 6.0 ) * CURV_r

                FLUX_DEN_x[k, 1, j] =  u_bnd[k, 1, j] * q_star - D_hor * ( GRAD_bnd_x[k, 1, j] - uΔt / 2.0 * CURV_r )
            end
        end
    end


    # y
    for i=1:Nx, j=2:Ny
        for k=1:Nz[i, j]
            if noflux_y_mask3[k, i, j] == 0.0
                FLUX_DEN_y[k, i, j] = 0.0
            else

                CURV_r = ( v_bnd[k, i, j] >= 0.0 ) ? CURV_y[k, i, j-1] : CURV_y[k, i, j]
                vΔt    = v_bnd[k, i, j] * Δt
                q_star = (qs[k, i, j-1] + qs[k, i, j]) / 2.0 - vΔt / 2.0 * GRAD_bnd_y[k, i, j] + ( D_hor * Δt / 2.0 - gi.dy_s[i, j]^2.0/6.0 + vΔt^2.0 / 6.0 ) * CURV_r

                FLUX_DEN_y[k, i, j] = v_bnd[k, i, j] * q_star - D_hor * ( GRAD_bnd_y[k, i, j] - vΔt / 2.0 * CURV_r )

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

        #local q_star
        for k=2:_Nz

            CURV_r = ( w_bnd[k, i, j] >= 0.0 ) ? CURV_z[k, i, j] : CURV_z[k-1, i, j]
            wΔt    = w_bnd[k, i, j] * Δt
            q_star = (qs[k, i, j] + qs[k-1, i, j]) / 2.0 - wΔt / 2.0 * GRAD_bnd_z[k, i, j] + ( D_ver * Δt / 2.0 - Δzs[k-1, i, j]^2.0/6.0 + wΔt^2.0 / 6.0 ) * CURV_r

            FLUX_DEN_z[k, i, j] = w_bnd[k, i, j] * q_star - D_ver * ( GRAD_bnd_z[k, i, j] - wΔt / 2.0 * CURV_r )

        end

        FLUX_DEN_z[_Nz+1, i, j] = FLUX_bot[i, j] = FLUX_DEN_z[_Nz, i, j]

        #println("(i, j) = ", (i, j), "; q_star = ", q_star, "; w_bnd = ", w_bnd[_Nz+1, i, j])

    end

end



function calGRAD_CURV!(;
    gi         :: PolelikeCoordinate.GridInfo,
    Nx         :: Integer,
    Ny         :: Integer,
    Nz         :: AbstractArray{Int64, 2},
    qs         :: AbstractArray{Float64, 3},     # ( Nz_bone   ,  Nx  , Ny   )
    GRAD_bnd_x :: AbstractArray{Float64, 3},     # ( Nz_bone   ,  Nx  , Ny   )
    GRAD_bnd_y :: AbstractArray{Float64, 3},     # ( Nz_bone   ,  Nx  , Ny+1 )
    GRAD_bnd_z :: AbstractArray{Float64, 3},     # ( Nz_bone+1 ,  Nx  , Ny   )
    CURV_x     :: AbstractArray{Float64, 3},     # ( Nz_bone   ,  Nx  , Ny   )
    CURV_y     :: AbstractArray{Float64, 3},     # ( Nz_bone   ,  Nx  , Ny   )
    CURV_z     :: AbstractArray{Float64, 3},     # ( Nz_bone   ,  Nx  , Ny   )
    mask3      :: AbstractArray{Float64, 3},     # ( Nz_bone   ,  Nx  , Ny   )
    noflux_x_mask3 :: AbstractArray{Float64, 3}, # ( Nz_bone   ,  Nx  , Ny   )
    noflux_y_mask3 :: AbstractArray{Float64, 3}, # ( Nz_bone   ,  Nx  , Ny+1 )
    Δzs        :: AbstractArray{Float64, 3},     # ( Nz_bone-1 ,  Nx  , Ny   )
    hs         :: AbstractArray{Float64, 3},     # ( Nz_bone   ,  Nx  , Ny   )
)

    # x 
    for i=2:Nx, j=1:Ny 
        for k=1:Nz[i, j]
            GRAD_bnd_x[k, i, j] = (
                ( noflux_x_mask3[k, i, j] == 0.0 )  
                ? 0.0 : ( qs[k, i, j] - qs[k, i-1, j] ) / gi.dx_w[i, j] 
            )
        end
    end

    # x - periodic
    for j=1:Ny
        for k=1:Nz[1, j]
            GRAD_bnd_x[k, 1, j] = (
                ( noflux_x_mask3[k, 1, j] == 0.0 )  
                ? 0.0 : ( qs[k, 1, j] - qs[k, Nx, j] ) / gi.dx_w[1, j]
            )
        end
    end

    # y
    for i=1:Nx, j=2:Ny
        for k=1:Nz[i, j]
            GRAD_bnd_y[k, i, j] = (
                ( noflux_y_mask3[k, i, j] == 0.0 )  
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

        GRAD_bnd_z[_Nz+1, i, j] = GRAD_bnd_z[_Nz, i, j]

    end

    # CURV
    for i=1:Nx, j=1:Ny
        for k=1:Nz[i, j]
            CURV_x[k, i, j] = ( GRAD_bnd_x[k, i+1, j  ] - GRAD_bnd_x[k  , i, j] ) / gi.dx_c[i, j]
            CURV_y[k, i, j] = ( GRAD_bnd_y[k, i  , j+1] - GRAD_bnd_y[k  , i, j] ) / gi.dy_c[i, j]
            CURV_z[k, i, j] = ( GRAD_bnd_z[k, i  , j  ] - GRAD_bnd_z[k+1, i, j] ) / hs[k, i, j]
#            if (k, i, j) == (4, 47, 87)
#               #println("[3,47,87] CURV_z=", CURV_z[1:6, i, j], ", hs=", hs[1:6, i, j], "; Δzs: ", Δzs[1:6, i, j])
#               #println("[3,47,87] GRAD_bnd_z: ", GRAD_bnd_z[1:6, i, j])
#            end
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


function calDIV!(;
    ASUM     ,
    u_bnd     :: AbstractArray{Float64, 3},   # ( Nz_bone  , Nx  , Ny   )
    v_bnd     :: AbstractArray{Float64, 3},   # ( Nz_bone  , Nx  , Ny+1 )
    div       :: AbstractArray{Float64, 3},   # ( Nz_bone  , Nx  , Ny   )
    workspace :: AbstractArray{Float64, 3},   # ( Nz_bone  , Nx  , Ny   )
)

    
    mul!(view(div, :),       ASUM.mtx_DIV_X, view(u_bnd, :))
    mul!(view(workspace, :), ASUM.mtx_DIV_Y, view(v_bnd, :))
    div .+= workspace

end

#=
function calDIV!(;
    gi       :: PolelikeCoordinate.GridInfo,
    Nx       :: Integer,
    Ny       :: Integer,
    Nz       :: AbstractArray{Int64, 2},
    u_bnd    :: AbstractArray{Float64, 3},   # ( Nz_bone  , Nx  , Ny   )
    v_bnd    :: AbstractArray{Float64, 3},   # ( Nz_bone  , Nx  , Ny+1 )
    div      :: AbstractArray{Float64, 3},   # ( Nz_bone  , Nx  , Ny   )
    mask3    :: AbstractArray{Float64, 3},   # ( Nz_bone  , Nx  , Ny   )
)

    
    mul!(view(div, :),           ASUM.mtx_DIV_X, view(u_bnd, :))
    mul!(view(workspaces[1], :), ASUM.mtx_DIV_Y, view(v_bnd, :))
    div .+= workspaces[1]

    #=

#    local tmp = tmp_σ = 0.0
    for i=1:Nx, j=1:Ny

        for k=1:Nz[i, j]

            if mask3[k, i, j] == 0.0
                break
            end
            
            div[k, i, j] =  (  
                u_bnd[k, i+1, j  ]  * gi.DY[i+1, j  ]
              - u_bnd[k, i,   j  ]  * gi.DY[i  , j  ]
              + v_bnd[k, i,   j+1]  * gi.DX[i  , j+1]
              - v_bnd[k, i,   j  ]  * gi.DX[i  , j  ]
            ) / gi.dσ[i, j]

        end

    end

    =#
end
=#


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

        w_bnd[1, i, j] = 0.0

        for k=1:Nz[i, j]

            if mask3[k, i, j] == 0.0
                break
            end
            
            w_bnd[k+1, i, j] = w_bnd[k, i, j] + div[k, i, j] * hs[k, i, j]
        end

#        tmp   += w_bnd[Nz[i, j]+1, i, j] * gi.dσ[i, j]
#        tmp_σ += gi.dσ[i, j]
    end

#   #println("tmp: ", tmp, "; tmp_σ: ", tmp_σ, "; Average w: ", tmp/tmp_σ)

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
    noflux_x_mask3 :: AbstractArray{Float64, 3}, # (Nz_bone, Nx+1, Ny)
    noflux_y_mask3 :: AbstractArray{Float64, 3}, # (Nz_bone, Nx, Ny+1)
)

    # x
    for i=2:Nx, j=1:Ny
        for k=1:Nz[i, j]
            if noflux_x_mask3[k, i, j] == 0.0
                u_bnd[k, i, j] = 0.0
            else
                u_bnd[k, i, j] = u[k, i-1, j] * (1.0 - weight_e[i, j]) + u[k, i, j] * weight_e[i, j]
                #u_bnd[k, i, j] = (u[k, i-1, j] + u[k, i, j]) / 2.0
            end
        end
    end
    
    # x - periodic
    for j=1:Ny
        for k=1:Nz[1, j]
            if noflux_x_mask3[k, 1, j] == 0.0
                u_bnd[k, 1, j] = u_bnd[k, Nx+1, j] = 0.0
            else
                u_bnd[k, 1, j] = u_bnd[k, Nx+1, j] = u[k, Nx, j] * (1.0 - weight_e[1, j]) + u[k, 1, j] * weight_e[1, j]
                #u_bnd[k, 1, j] = u_bnd[k, Nx+1, j] = (u[k, Nx, j] + u[k, 1, j]) / 2.0
            end
        end
    end

    # y
    for i=1:Nx, j=2:Ny
        for k=1:Nz[i, j]
            if noflux_y_mask3[k, i, j] == 0.0
                v_bnd[k, i, j] = 0.0
            else
                v_bnd[k, i, j] = v[k, i, j-1] * (1.0 - weight_n[i, j]) + v[k, i, j] * weight_n[i, j]
                #v_bnd[k, i, j] = (v[k, i, j-1] + v[k, i, j]) / 2.0
            end
        end
    end

end
