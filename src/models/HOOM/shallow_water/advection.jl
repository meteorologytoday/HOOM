
function calDiffAdv_QUICK_SpeedUp!(
    model       :: Model,
    Δt          :: Float64,
)
   
    tcr_adv = model.tcr_adv
    env     = model.env
    state   = model.state
    ASUM    = tcr_adv.ASUM
    workspaces = tcr_adv.workspaces

    for x=1:env.NX
 
        X            = view(model.state.X,        :, :, :, x)
        XFLUX_bot    = view(tcr_adv.XFLUX_bot,       :, :, x)
        XFLUX_CONV   = view(tcr_adv.XFLUX_CONV,   :, :, :, x)
        XFLUX_CONV_h = view(tcr_adv.XFLUX_CONV_h, :, :, :, x)
        XFLUX_DEN_x  = view(tcr_adv.XFLUX_DEN_x,  :, :, :, x)
        XFLUX_DEN_y  = view(tcr_adv.XFLUX_DEN_y,  :, :, :, x)
        XFLUX_DEN_z  = view(tcr_adv.XFLUX_DEN_z,  :, :, :, x)

        let
            mul!(view(tcr_adv.GRAD_bnd_x, :), ASUM.mtx_GRAD_X, view(X, :))
            mul!(view(tcr_adv.GRAD_bnd_y, :), ASUM.mtx_GRAD_Y, view(X, :))
            mul!(view(tcr_adv.GRAD_bnd_z, :), ASUM.mtx_GRAD_Z, view(X, :))
     
            mul!(view(tcr_adv.CURV_x, :), ASUM.mtx_CURV_X, view(tcr_adv.GRAD_bnd_x, :))
            mul!(view(tcr_adv.CURV_y, :), ASUM.mtx_CURV_Y, view(tcr_adv.GRAD_bnd_y, :))
            mul!(view(tcr_adv.CURV_z, :), ASUM.mtx_CURV_Z, view(tcr_adv.GRAD_bnd_z, :))
     
        end

   #println("FLUXDEN")
        calFluxDensity!(
            gi         = env.gi,
            Nx         = env.Nx,
            Ny         = env.Ny,
            Nz         = env.Nz_av_f,
            FLUX_bot   = XFLUX_bot,
            qs         = X,
            GRAD_bnd_x = tcr_adv.GRAD_bnd_x,
            GRAD_bnd_y = tcr_adv.GRAD_bnd_y,
            GRAD_bnd_z = tcr_adv.GRAD_bnd_z,
            CURV_x     = tcr_adv.CURV_x,
            CURV_y     = tcr_adv.CURV_y,
            CURV_z     = tcr_adv.CURV_z,
            FLUX_DEN_x = XFLUX_DEN_x,
            FLUX_DEN_y = XFLUX_DEN_y,
            FLUX_DEN_z = XFLUX_DEN_z,
            u_bnd      = state.u_f,
            v_bnd      = state.v_f,
            w_bnd      = state.w_f,
            mask3          = env.mask3_f,
            noflux_x_mask3 = env.noflux_x_mask3_f,
            noflux_y_mask3 = env.noflux_y_mask3_f,
            Δzs        = env.Δz_f,
            D_hor      = env.Dh[x],
            D_ver      = env.Dv[x],
            Δt         = Δt,
        )


   #println("TOTAL CHANGE")
        let
            mul!(view(XFLUX_CONV_h, :), ASUM.mtx_DIV_X, view(XFLUX_DEN_x, :))
            mul!(view(workspaces[2], :),   ASUM.mtx_DIV_Y, view(XFLUX_DEN_y, :))
            mul!(view(workspaces[3], :),   ASUM.mtx_DIV_Z, view(XFLUX_DEN_z, :))

            XFLUX_CONV_h .+= workspaces[2]
            XFLUX_CONV_h .*= -1.0
            XFLUX_CONV .= XFLUX_CONV_h 
            XFLUX_CONV .-= workspaces[3]

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
    gi         :: DisplacedPoleCoordinate.GridInfo,
    Nx         :: Integer,
    Ny         :: Integer,
    Nz         :: AbstractArray{Int64, 2},
    FLUX_bot     :: AbstractArray{Float64, 2},     # ( Nx  , Ny   )
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
    noflux_x_mask3 :: AbstractArray{Float64, 3}, # ( Nz_bone   ,  Nx+1, Ny   )
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
                FLUX_DEN_x[k, 1, j] = FLUX_DEN_x[k, Nx+1, j] = 0.0
            else

                CURV_r = ( u_bnd[k, 1, j] >= 0.0 ) ? CURV_x[k, Nx, j] : CURV_x[k, 1, j]
                uΔt    = u_bnd[k, 1, j] * Δt
                q_star = (qs[k, Nx, j] + qs[k, 1, j]) / 2.0 - uΔt / 2.0 * GRAD_bnd_x[k, 1, j] + ( D_hor * Δt / 2.0 - gi.dx_w[1, j]^2.0/6.0 + uΔt^2.0 / 6.0 ) * CURV_r

                FLUX_DEN_x[k, 1, j] = FLUX_DEN_x[k, Nx+1, j] = u_bnd[k, 1, j] * q_star - D_hor * ( GRAD_bnd_x[k, 1, j] - uΔt / 2.0 * CURV_r )
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
    noflux_x_mask3 :: AbstractArray{Float64, 3}, # ( Nz_bone   ,  Nx+1, Ny   )
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
            GRAD_bnd_x[k, 1, j] = GRAD_bnd_x[k, Nx+1, j] = (
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
    gi       :: DisplacedPoleCoordinate.GridInfo,
    Nx       :: Integer,
    Ny       :: Integer,
    Nz       :: AbstractArray{Int64, 2},
    u_bnd    :: AbstractArray{Float64, 3},   # ( Nz_bone  , Nx+1, Ny   )
    v_bnd    :: AbstractArray{Float64, 3},   # ( Nz_bone  , Nx  , Ny+1 )
    div      :: AbstractArray{Float64, 3},   # ( Nz_bone  , Nx  , Ny   )
    mask3    :: AbstractArray{Float64, 3},   # ( Nz_bone  , Nx  , Ny   )
)

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


end



function calVerVelBnd!(;
    gi       :: DisplacedPoleCoordinate.GridInfo,
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
