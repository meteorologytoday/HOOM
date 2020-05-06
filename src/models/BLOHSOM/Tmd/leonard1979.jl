function calDIV!(;
    ASUM      :: AdvectionSpeedUpMatrix,
    u_U       :: AbstractArray{Float64, 3},   # ( Nz_bone  , Nx  , Ny   )
    v_V       :: AbstractArray{Float64, 3},   # ( Nz_bone  , Nx  , Ny+1 )
    div       :: AbstractArray{Float64, 3},   # ( Nz_bone  , Nx  , Ny   )
    workspace :: AbstractArray{Float64, 3},   # ( Nz_bone  , Nx  , Ny   )
)

    
    mul_autoflat!(div,       ASUM.T_DIVx_U, u_U)
    mul_autoflat!(workspace, ASUM.T_DIVy_V, v_V)
    @. div += workspace

end


function calDiffAdv_QUICKEST_SpeedUp!(
    m       :: TmdModel,
    Δt          :: Float64,
)

    @fast_extract m
    
    ASUM    = co.ASUM
    wksp    = co.wksp

    GRADx_U = getSpace!(wksp, :U)
    GRADy_V = getSpace!(wksp, :V)
    GRADz_W = getSpace!(wksp, :W)
        
    CURVx_T = getSpace!(wksp, :T)
    CURVy_T = getSpace!(wksp, :T)
    CURVz_T = getSpace!(wksp, :T)
        
    tmp1 = getSpace!(wksp, :T)
    tmp2 = getSpace!(wksp, :T)
    tmp3 = getSpace!(wksp, :T)

    for x=1:ev.NX
 
        X            = view(st.X,            :, :, :, x)
        dΔXdt        = view(st.dΔXdt,           :, :, x)
        XFLUX_bot    = view(co.XFLUX_bot,       :, :, x)
        XFLUX_CONV   = view(co.XFLUX_CONV,   :, :, :, x)
        XFLUX_CONV_h = view(co.XFLUX_CONV_h, :, :, :, x)
        XFLUX_DEN_x  = view(co.XFLUX_DEN_x,  :, :, :, x)
        XFLUX_DEN_y  = view(co.XFLUX_DEN_y,  :, :, :, x)
        XFLUX_DEN_z  = view(co.XFLUX_DEN_z,  :, :, :, x)

        mul_autoflat!(GRADx_U, ASUM.U_∂x_T, X)
        mul_autoflat!(GRADy_V, ASUM.V_∂y_T, X)
        mul_autoflat!(GRADz_W, ASUM.W_∂z_T, X)

        mul_autoflat!(CURVx_T, ASUM.T_∂x_U, GRADx_U)
        mul_autoflat!(CURVy_T, ASUM.T_∂y_V, GRADy_V)
        mul_autoflat!(CURVz_T, ASUM.T_∂z_W, GRADz_W)

        Kh = ev.Kh_X[x]
        Kv = ev.Kv_X[x]

        # x direction
        calFluxDensity_abstract!(;
            X_T                 = X,
            u_bnd               = fr.u_U,
            GRAD_bnd            = GRADx_U,
            CURV_T              = CURVx_T,
            FLUX_DEN_bnd        = XFLUX_DEN_x,
            K                   = Kh,
            Δt                  = Δt,
            Δx_bnd              = ASUM.Δx_U,
            
            op_bnd_∂x_T         = ASUM.U_∂x_T,
            op_T_∂x_bnd         = ASUM.T_∂x_U,
            op_bnd_pos_dir_T    = ASUM.op.U_E_T,
            op_bnd_neg_dir_T    = ASUM.op.U_W_T,
            op_bnd_interp_T     = ASUM.U_interp_T,
            op_filter_bnd       = ASUM.filter_U,

            pos_u_mask_bnd      = getSpace!(co.wksp, :U),
            CURV_on_pos_dir_bnd = getSpace!(co.wksp, :U),
            CURV_on_neg_dir_bnd = getSpace!(co.wksp, :U),
            CURV_r_bnd          = getSpace!(co.wksp, :U),
            uΔt_bnd             = getSpace!(co.wksp, :U),
            X_star_bnd          = getSpace!(co.wksp, :U),
        )

        # y direction
        calFluxDensity_abstract!(;
            X_T                 = X,
            u_bnd               = fr.v_V,
            GRAD_bnd            = GRADy_V,
            CURV_T              = CURVy_T,
            FLUX_DEN_bnd        = XFLUX_DEN_y,
            K                   = Kh,
            Δt                  = Δt,
            Δx_bnd              = ASUM.Δy_V,
            
            op_bnd_∂x_T         = ASUM.V_∂y_T,
            op_T_∂x_bnd         = ASUM.T_∂y_V,
            op_bnd_pos_dir_T    = ASUM.op.V_N_T,
            op_bnd_neg_dir_T    = ASUM.op.V_S_T,
            op_bnd_interp_T     = ASUM.V_interp_T,
            op_filter_bnd       = ASUM.filter_V,

            pos_u_mask_bnd      = getSpace!(co.wksp, :V),
            CURV_on_pos_dir_bnd = getSpace!(co.wksp, :V),
            CURV_on_neg_dir_bnd = getSpace!(co.wksp, :V),
            CURV_r_bnd          = getSpace!(co.wksp, :V),
            uΔt_bnd             = getSpace!(co.wksp, :V),
            X_star_bnd          = getSpace!(co.wksp, :V),
        )

        # z direction
        calFluxDensity_abstract!(;
            X_T                 = X,
            u_bnd               = fr.w_W,
            GRAD_bnd            = GRADz_W,
            CURV_T              = CURVz_T,
            FLUX_DEN_bnd        = XFLUX_DEN_z,
            K                   = Kv,
            Δt                  = Δt,
            Δx_bnd              = ASUM.Δz_W,
            
            op_bnd_∂x_T         = ASUM.W_∂z_T,
            op_T_∂x_bnd         = ASUM.T_∂z_W,
            op_bnd_pos_dir_T    = ASUM.op.W_UP_T,
            op_bnd_neg_dir_T    = ASUM.op.W_DN_T,
            op_bnd_interp_T     = ASUM.W_interp_T,
            op_filter_bnd       = ASUM.filter_W,

            pos_u_mask_bnd      = getSpace!(co.wksp, :W),
            CURV_on_pos_dir_bnd = getSpace!(co.wksp, :W),
            CURV_on_neg_dir_bnd = getSpace!(co.wksp, :W),
            CURV_r_bnd          = getSpace!(co.wksp, :W),
            uΔt_bnd             = getSpace!(co.wksp, :W),
            X_star_bnd          = getSpace!(co.wksp, :W),
        )
 
        mul_autoflat!(XFLUX_CONV_h, ASUM.T_DIVx_U, XFLUX_DEN_x)
        mul_autoflat!(tmp1,   ASUM.T_DIVy_V, XFLUX_DEN_y)
        mul_autoflat!(tmp2,   ASUM.T_DIVz_W, XFLUX_DEN_z)


        @. XFLUX_CONV_h = -1.0 * (XFLUX_CONV_h + tmp1)
        @. XFLUX_CONV = XFLUX_CONV_h - tmp2

        
        calMixedLayer_dΔqdt!(
            Nx          = ev.Nx,
            Ny          = ev.Ny,
            Nz          = ev.Nz_av,
            FLUX_CONV_h = XFLUX_CONV_h,
            FLUX_DEN_z  = XFLUX_DEN_z,
            dΔqdt       = dΔXdt,
            mask        = ev.mask2,
            FLDO        = st.FLDO,
            h_ML        = st.h_ML,
            hs          = co.Δz_T,
            zs          = ev.z_bnd_av,
        )

    end

end

function calFluxDensity_abstract!(;
    X_T                       :: AbstractArray{Float64, 3},
    u_bnd                     :: AbstractArray{Float64, 3},
    GRAD_bnd                  :: AbstractArray{Float64, 3}, 
    CURV_T                    :: AbstractArray{Float64, 3},
    FLUX_DEN_bnd              :: AbstractArray{Float64, 3},
    K                         :: Float64,
    Δt                        :: Float64,
    Δx_bnd                    :: AbstractArray{Float64, 3},

    op_bnd_∂x_T               :: AbstractArray{Float64, 2},
    op_T_∂x_bnd               :: AbstractArray{Float64, 2},
    op_bnd_pos_dir_T          :: AbstractArray{Float64, 2},
    op_bnd_neg_dir_T          :: AbstractArray{Float64, 2},
    op_bnd_interp_T           :: AbstractArray{Float64, 2},
    op_filter_bnd             :: AbstractArray{Float64, 2},

    pos_u_mask_bnd            :: AbstractArray{Float64, 3},
    CURV_on_pos_dir_bnd       :: AbstractArray{Float64, 3},
    CURV_on_neg_dir_bnd       :: AbstractArray{Float64, 3},
    CURV_r_bnd                :: AbstractArray{Float64, 3},
    uΔt_bnd                   :: AbstractArray{Float64, 3},
    X_star_bnd                :: AbstractArray{Float64, 3},
)

    mul_autoflat!(GRAD_bnd, op_bnd_∂x_T, X_T)
    mul_autoflat!(CURV_T, op_T_∂x_bnd, GRAD_bnd)

    @. pos_u_mask_bnd = u_bnd >= 0
    mul_autoflat!(CURV_on_pos_dir_bnd, op_bnd_neg_dir_T, CURV_T)
    mul_autoflat!(CURV_on_neg_dir_bnd, op_bnd_pos_dir_T, CURV_T)
    @. CURV_r_bnd = CURV_on_neg_dir_bnd * pos_u_mask_bnd + CURV_on_pos_dir_bnd * (1-pos_u_mask_bnd)

    @. uΔt_bnd    = u_bnd * Δt
    mul_autoflat!(X_star_bnd, op_bnd_interp_T, X_T)
    @. X_star_bnd = X_star_bnd - uΔt_bnd / 2.0 * GRAD_bnd + ( K * Δt / 2.0 - Δx_bnd^2.0/6.0 + uΔt_bnd^2.0 / 6.0 ) * CURV_r_bnd

    # reuse
    empty = uΔt_bnd
    @. empty = u_bnd * X_star_bnd - K * ( GRAD_bnd - uΔt_bnd / 2.0 * CURV_r_bnd )

    mul_autoflat!(FLUX_DEN_bnd, op_filter_bnd, empty)

    

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

