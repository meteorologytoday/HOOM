function decompose!(;
    total   :: AbstractArray{Float64, 3},
    mean    :: AbstractArray{Float64, 2},
    anomaly :: AbstractArray{Float64, 3},
    mtx     :: AvgSpeedUpMatrix,
    mask2   :: AbstractArray{Float64, 2},
)

    Nx, Ny = size(mask2)

    mul2!(view(mean, :), mtx.S_avg_C, view(total, :))

    for i=1:Nx, j=1:Ny

        if mask2[i, j] == 0
            continue
        end

        total[:, i, j] .-= mean[i, j]

    end
end

function advectDynamic!(
    model   :: Model,
    Δt      :: Float64,
)
    calG!(model, dyn_adv.G_idx[:now])
end

function calG!(
    model   :: Model,
    G_idx   :: Int64,
)

    # 1. derive barotropic and baroclinic flow
    # 2. derive G of each terms
    #
    #    a)
    #
    state = model.state
    dyn_adv = model.dyn_adv
    ASUM = dyn_adv.ASUM
    AVGM = dyn_adv.AVGM
    env = model.env

    # cal b_f from T_f, S_f
    for i=1:env.Nx, j=1:env.Ny

        if env.mask2[i, j] == 0
            continue
        end

        for k=1:env.Nz_av_f[i, j]
            state.b_f[k, i, j] = TS2b(state.T[k, i, j], state.S[k, i, j])
        end
    end

    # cal b_c from b_f
    mul2!(state.b_c, AVGM.C_avg_F, state.b_f) 
    
    # cal barotropic and baroclinic components
    decompose!(total=state.u_c, mean=state.U, anomaly=state.u, mtx=AVGM)
    decompose!(total=state.v_c, mean=state.V, anomaly=state.v, mtx=AVGM)
    decompose!(total=state.b_c, mean=state.B, anomaly=state.b, mtx=AVGM)
  
    # cal τx_acc, τy_acc
    mul2!(dyn_adv.τx_acc, ASUM.U_interp_T, state.τx)
    mul2!(dyn_adv.τy_acc, ASUM.V_interp_T, state.τy)
    dyn_adv.τx_acc /= (ρ_sw * env.H_c[1])
    dyn_adv.τy_acc /= (ρ_sw * env.H_c[1])
 
    # cal∇Φ
    mul2!(dyn_adv.∂Φ∂x, ASUM.U_∂x_T, state.Φ)
    mul2!(dyn_adv.∂Φ∂y, ASUM.U_∂y_T, state.Φ)
 
    # cal ∇b
    mul2!(dyn_adv.∂b∂x, ASUM.U_∂x_T, state.b)
    mul2!(dyn_adv.∂b∂y, ASUM.V_∂y_T, state.b)

    # cal Coriolis force
    mul2!(dyn_adv.fu, ASUM.V_f_U, state.u)   # fu on V grid
    mul2!(dyn_adv.fv, ASUM.U_f_V, state.v)   # fv on V grid
    
    # ===== [ BEGIN cal (v⋅∇)v ] =====
    mul2!(dyn_adv.∂u∂x, ASUM.U_∂x_U, state.u)
    mul2!(dyn_adv.∂u∂y, ASUM.U_∂y_U, state.u)
    mul2!(dyn_adv.∂v∂x, ASUM.V_∂x_V, state.v)
    mul2!(dyn_adv.∂v∂y, ASUM.V_∂y_V, state.v)
    
    # interpolate first then multiply by ∇v
    
    # On U grid
    mul2!(dyn_adv.v∂u∂y, ASUM.U_interp_V, state.v_c)  # store interpolated v into v∂u∂y
    for k=1:env.Nz_c, i=1:env.Nx, j=1:env.Ny
        dyn_adv.u∂u∂x[k, i, j]  = state.u_c[k, i, j] * dyn_adv.∂u∂x[k, i, j]
        dyn_adv.v∂u∂y[k, i, j] *=                      dyn_adv.∂u∂y[k, i, j]
    end

    # On V grid
    mul2!(dyn_adv.u∂v∂x, ASUM.V_interp_U, state.u_c)  # store interpolated u into u∂v∂x
    for k=1:env.Nz_c, i=1:env.Nx, j=1:env.Ny+1
        dyn_adv.u∂v∂x[k, i, j] *=                      dyn_adv.∂v∂x[k, i, j]
        dyn_adv.v∂v∂y[k, i, j]  = state.v_c[k, i,j]  * dyn_adv.∂v∂y[k, i, j]
    end
    # ===== [ END cal (v⋅∇)v ] =====

    # cal G_bc

    G_bc_u = view(dyn_adv.G_bc_u, :, :, :, G_idx)
    G_bc_v = view(dyn_adv.G_bc_v, :, :, :, G_idx)
    G_bt_u = view(dyn_adv.G_bt_u, :, :, :, G_idx)
    G_bt_v = view(dyn_adv.G_bt_v, :, :, :, G_idx)

    G_bc_u .= 0 
    G_bc_u -= dyn_adv.u∂u∂x
    G_bc_u -= dyn_adv.v∂u∂y
    G_bc_u += dyn_adv.τx_acc
    G_bc_u -= dyn_adv.∂Φ∂x
    G_bc_u += dyn_adv.∂b∂x
    G_bc_u += dyn_adv.fv
 
    G_bc_v .= 0 
    G_bc_v -= dyn_adv.u∂v∂x
    G_bc_v -= dyn_adv.v∂v∂y
    G_bc_v += dyn_adv.τy_acc
    G_bc_v -= dyn_adv.∂Φ∂y
    G_bc_v += dyn_adv.∂b∂x
    G_bc_v -= dyn_adv.fu
    
end



