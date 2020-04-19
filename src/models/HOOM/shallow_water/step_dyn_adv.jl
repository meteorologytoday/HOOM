function decompose!(;
    total   :: AbstractArray{Float64, 3},
    mean    :: AbstractArray{Float64, 2},
    anomaly :: AbstractArray{Float64, 3},
    mask    :: AbstractArray{Float64, 2},
    va       :: VerticalAverager,
)

    Nx, Ny = size(mask)

    
    calAverage_c2s!(
        va,
        total,
        mean;
        mask=mask,
    )

    for i=1:Nx, j=1:Ny

        if mask[i, j] == 0
            continue
        end

        total[:, i, j] .-= mean[i, j]

    end
end

function advectDynamic!(
    model   :: DynModel,
)
    env   = model.env
    state = model.state
    core  = model.core

    reset!(core.wksp)

    #println("calAuxV!")
    calAuxV!(model)
    #println("calAuxΦ!")
    calAuxΦ!(model)
    #println("solveΦ!")
    solveΦ!(model)
    #println("updateV!")
    updateV!(model)

#    projectVertical!(va, )
  
end

function calAuxV!(
    model   :: DynModel,
)

    # 1. derive barotropic and baroclinic flow
    # 2. derive G of each terms
    #
    #    a)
    #
    state = model.state
    core = model.core
    c_ops = core.c_ops
    s_ops = core.s_ops
    va = core.va
    env = model.env
    wksp = core.wksp

    # cal b_f from T_f, S_f
    for i=1:env.Nx, j=1:env.Ny

        if env.mask[i, j] == 0
            continue
        end

        for k=1:env.Nz_f
            state.b_f[k, i, j] = TS2b(state.T[k, i, j], state.S[k, i, j])
        end
    end

    # cal b_c from b_f
    calAverage_f2c!(va, state.b_f, state.b_c; mask=env.mask)
    
    # cal barotropic and baroclinic components
    decompose!(total=state.u_c, mean=state.U, anomaly=state.u, va=va, mask=env.mask)
    decompose!(total=state.v_c, mean=state.V, anomaly=state.v, va=va, mask=env.mask)
    decompose!(total=state.b_c, mean=state.B, anomaly=state.b, va=va, mask=env.mask)

  
    # cal τx_acc, τy_acc
    τx_acc = getSpace!(wksp, :sU)
    τy_acc = getSpace!(wksp, :sV)

    mul2!(τx_acc, s_ops.U_interp_T, state.τx)
    mul2!(τy_acc, s_ops.V_interp_T, state.τy)
    τx_acc /= (ρ_sw * env.H_c[1])
    τy_acc /= (ρ_sw * env.H_c[1])
 
    # cal ∇b
    ∂b∂x   = getSpace!(wksp, :cU)
    ∂b∂y   = getSpace!(wksp, :cV)
    mul2!(∂b∂x, c_ops.U_∂x_T, state.b)
    mul2!(∂b∂y, c_ops.V_∂y_T, state.b)

    # cal Coriolis force
    fu   = getSpace!(wksp, :cV)
    fv   = getSpace!(wksp, :cU)
    mul2!(fu, c_ops.V_f_U, state.u)   # fu on V grid
    mul2!(fv, c_ops.U_f_V, state.v)   # fv on U grid
    
    # ===== [ BEGIN cal (v⋅∇)v ] =====
    ∂u∂x = getSpace!(wksp, :cU)
    ∂u∂y = getSpace!(wksp, :cU)
    ∂v∂x = getSpace!(wksp, :cV)
    ∂v∂y = getSpace!(wksp, :cV)
    mul2!(∂u∂x, c_ops.U_∂x_U, state.u)
    mul2!(∂u∂y, c_ops.U_∂y_U, state.u)
    mul2!(∂v∂x, c_ops.V_∂x_V, state.v)
    mul2!(∂v∂y, c_ops.V_∂y_V, state.v)
    
    # interpolate first then multiply by ∇v
    
    # On U grid
    u∂u∂x = getSpace!(wksp, :cU)
    v∂u∂y = getSpace!(wksp, :cU)
    mul2!(v∂u∂y, c_ops.U_interp_V, state.v_c)  # store interpolated v into v∂u∂y
    for k=1:env.Nz_c, i=1:env.Nx, j=1:env.Ny
        u∂u∂x[k, i, j]  = state.u_c[k, i, j] * ∂u∂x[k, i, j]
        v∂u∂y[k, i, j] *=                      ∂u∂y[k, i, j]
    end

    # On V grid
    u∂v∂x = getSpace!(wksp, :cV)
    v∂v∂y = getSpace!(wksp, :cV)
    mul2!(u∂v∂x, c_ops.V_interp_U, state.u_c)  # store interpolated u into u∂v∂x
    for k=1:env.Nz_c, i=1:env.Nx, j=1:env.Ny+1
        u∂v∂x[k, i, j] *=                      ∂v∂x[k, i, j]
        v∂v∂y[k, i, j]  = state.v_c[k, i,j]  * ∂v∂y[k, i, j]
    end
    # ===== [ END cal (v⋅∇)v ] =====

    # cal G
    G_idx = core.G_idx
    Δt0 = G_idx[:now]
    Δt1 = G_idx[:one_Δt_ago]
    Δt2 = G_idx[:two_Δt_ago]


    G_u = view(state.G_u, :, :, :, Δt0)
    G_v = view(state.G_v, :, :, :, Δt0)

    G_u .= 0 
    G_u -= u∂u∂x
    G_u -= v∂u∂y
    G_u += ∂b∂x
    G_u += fv


 
    G_v .= 0 
    G_v -= u∂v∂x
    G_v -= v∂v∂y
    G_v += ∂b∂y
    G_v -= fu

    # surface
    G_u[1, :, :] += τx_acc
    G_v[1, :, :] += τy_acc

    # calculate auxiliary velocity
    Δt = env.Δt 
    G_u = state.G_u
    G_v = state.G_v
    for k=1:env.Nz_c, i=1:env.Nx, j=1:env.Ny
        core.u_aux[k, i, j] = state.u_c[k, i, j] + Δt *
           ABIII(
                G_u[k, i, j, Δt0],
                G_u[k, i, j, Δt1],
                G_u[k, i, j, Δt2],
        )
    end

    for k=1:env.Nz_c, i=1:env.Nx, j=1:env.Ny+1
        core.v_aux[k, i, j] = state.v_c[k, i, j] + Δt *
            ABIII(
                G_v[k, i, j, Δt0],
                G_v[k, i, j, Δt1],
                G_v[k, i, j, Δt2],
        )
       
    end

    G_idx[:now] = Δt2
    G_idx[:one_Δt_ago] = Δt0
    G_idx[:two_Δt_ago] = Δt1

end

function calAuxΦ!(
    model :: DynModel,
)
    # cal mean aux_v
    core = model.core
    env  = model.env
    state = model.state

    s_ops = core.s_ops
    va    = core.va
    wksp  = core.wksp
    mask  = env.mask
    Δt    = env.Δt

    
    Hu     = getSpace!(wksp, :sU)
    Hv     = getSpace!(wksp, :sV)
    DIV_Hu = getSpace!(wksp, :sT)
    DIV_Hv = getSpace!(wksp, :sT)

    calAverage_c2s!(va, core.u_aux, Hu; mask=mask)
    calAverage_c2s!(va, core.v_aux, Hv; mask=mask)
    
    Hu .*= env.H_total
    Hv .*= env.H_total

    mul2!(DIV_Hu, s_ops.T_DIVx_U,   Hu)
    mul2!(DIV_Hv, s_ops.T_DIVy_V,   Hv)
 
    for i=1:env.Nx, j=1:env.Ny
        core.Φ_aux[i, j] = state.Φ[i, j] - Δt * (
            DIV_Hu[i, j] + DIV_Hv[i, j]
        )
    end
   
end

function solveΦ!(
    model :: DynModel,
)
    env    = model.env
    core   = model.core
    state  = model.state

    rhs = getSpace!(core.wksp, :sT)
    solver = core.Φ_solver
    α = solver.α
    for i=1:env.Nx, j=1:env.Ny
        rhs[i, j] = - core.Φ_aux[i, j] * α
    end
 
    ldiv!(
        view(state.Φ, :),
        solver.tool_mtx.MoLap,
        view(rhs, :),
    )

end

function updateV!(
    model :: DynModel,
)
   
    core  = model.core
    env   = model.env
    state = model.state

    s_ops = core.s_ops
    wksp = core.wksp
    Δt   = model.env.Δt
    
    Δt∂Φ∂x = getSpace!(wksp, :sU)
    Δt∂Φ∂y = getSpace!(wksp, :sV)
     
    mul2!(Δt∂Φ∂x, s_ops.U_∂x_T, state.Φ)
    mul2!(Δt∂Φ∂y, s_ops.V_∂y_T, state.Φ)

    Δt∂Φ∂x .*= Δt
    Δt∂Φ∂y .*= Δt
 
    for k=1:env.Nz_c, i=1:env.Nx, j=1:env.Ny
        state.u_c[k, i, j] = core.u_aux[k, i, j] - Δt∂Φ∂x[i, j]
        state.v_c[k, i, j] = core.v_aux[k, i, j] - Δt∂Φ∂y[i, j]
    end
 
    projVertical_c2f!(core.va, state.u_c, state.u_f, mask=env.mask)
    projVertical_c2f!(core.va, state.v_c, state.v_f, mask=env.mask)
   
end

