function decompose!(;
    total   :: AbstractArray{Float64, 3},
    mean    :: AbstractArray{Float64, 2},
    anomaly :: AbstractArray{Float64, 3},
    va       :: VerticalAverager,
)

    Nx, Ny = size(mean)

    
    calAverage_c2s!(
        va,
        total,
        mean;
    )

    for i=1:Nx, j=1:Ny

        m = mean[i, j]
        for k=1:va.Nz_c
            anomaly[i, j, k] = total[i, j, k] - m
        end

    end
end

function advectDynamic!(
    model   :: DynModel,
)
    env   = model.env
    state = model.state
    core  = model.core



    #println("calAuxV!")
    calAuxV!(model)

    #println("calAuxΦ!")
    calAuxΦ!(model)

    #println("solveΦ!")
    solveΦ!(model)

    #println("updateV!")
    updateV!(model)

    #println(format("(u, v) = ({:.2f}, {:.2f})", state.u_c[1, 5, 5], state.v_c[1,5,5]))
end

function calAuxV!(
    model   :: DynModel,
)

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
            state.b_f[i, j, k] = TS2b(state.T[i, j, k], state.S[i, j, k])
        end
    end

    # cal b_c from b_f
    #@time let
    #println("cal avg f2c")
    #calAverage_f2c!(va, state.b_f, state.b_c)
    #end
    
    #println(format("u_c, U, u before {:.2f}, {:.2f}, {:.2f}", state.u_c[5,5,1], state.U[5,5], state.u[5,5,1]))
    # cal barotropic and baroclinic components
    decompose!(total=state.u_c, mean=state.U, anomaly=state.u, va=va)
    decompose!(total=state.v_c, mean=state.V, anomaly=state.v, va=va)
    decompose!(total=state.b_c, mean=state.B, anomaly=state.b, va=va)

    #println(format("u_c, U, u after {:.2f}, {:.2f}, {:.2f}", state.u_c[5,5,1], state.U[5,5], state.u[5,5,1]))
  
    # cal τx_acc, τy_acc
    τx_acc = getSpace!(wksp, :sU)
    τy_acc = getSpace!(wksp, :sV)

    mul2!(τx_acc, s_ops.U_interp_T, state.τx)
    mul2!(τy_acc, s_ops.V_interp_T, state.τy)
    τx_acc ./= (ρ_fw * env.H_c[1])
    τy_acc ./= (ρ_fw * env.H_c[1])
 
    # cal ∇b
    ∂b∂x   = getSpace!(wksp, :cU)
    ∂b∂y   = getSpace!(wksp, :cV)
    mul3!(∂b∂x, c_ops.U_∂x_T, state.b)
    mul3!(∂b∂y, c_ops.V_∂y_T, state.b)

    # cal Coriolis force
    fu   = getSpace!(wksp, :cV)
    fv   = getSpace!(wksp, :cU)
    mul3!(fu, c_ops.V_f_U, state.u_c)   # fu on V grid
    mul3!(fv, c_ops.U_f_V, state.v_c)   # fv on U grid
    #println("fu: ", fu[30, 29, 1])   

 
    # ===== [ BEGIN cal (v⋅∇)v ] =====
    ∂u∂x = getSpace!(wksp, :cU)
    ∂u∂y = getSpace!(wksp, :cU)
    ∂v∂x = getSpace!(wksp, :cV)
    ∂v∂y = getSpace!(wksp, :cV)
    mul3!(∂u∂x, c_ops.U_∂x_U, state.u)
    mul3!(∂u∂y, c_ops.U_∂y_U, state.u)
    mul3!(∂v∂x, c_ops.V_∂x_V, state.v)
    mul3!(∂v∂y, c_ops.V_∂y_V, state.v)
    
    # interpolate first then multiply by ∇v
    
    # On U grid
    u∂u∂x = getSpace!(wksp, :cU)
    v∂u∂y = getSpace!(wksp, :cU)
    mul3!(v∂u∂y, c_ops.U_interp_V, state.v_c)  # store interpolated v into v∂u∂y
    for i=1:env.Nx, j=1:env.Ny, k=1:env.Nz_c
        u∂u∂x[i, j, k]  = state.u[i, j, k] * ∂u∂x[i, j, k]
        v∂u∂y[i, j, k] *=                    ∂u∂y[i, j, k]
    end

    # On V grid
    u∂v∂x = getSpace!(wksp, :cV)
    v∂v∂y = getSpace!(wksp, :cV)
    mul3!(u∂v∂x, c_ops.V_interp_U, state.u_c)  # store interpolated u into u∂v∂x
    for i=1:env.Nx, j=1:env.Ny+1, k=1:env.Nz_c
        u∂v∂x[i, j, k] *=                     ∂v∂x[i, j, k]
        v∂v∂y[i, j, k]  = state.v[i, j, k]  * ∂v∂y[i, j, k]
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
#    G_u .-= u∂u∂x
#    G_u .-= v∂u∂y
#    G_u .+= ∂b∂x
    G_u .+= fv


 
    G_v .= 0 
#    G_v .-= u∂v∂x
#    G_v .-= v∂v∂y
#    G_v .+= ∂b∂y
    G_v .-= fu


    # surface
#    G_u[:, :, 1] .+= τx_acc
#    G_v[:, :, 1] .+= τy_acc
    #println("G_u")

    #=
    println("U")
    println(state.U)

    println("u first layer")
    println(state.u[1, :, :])

    println("dudx")
    println(∂u∂x[1, :, :])
    =#

    # calculate auxiliary velocity
    Δt = env.Δt 
    G_u = state.G_u
    G_v = state.G_v
    for i=1:env.Nx, j=1:env.Ny, k=1:env.Nz_c
        core.u_aux[i, j, k] = state.u_c[i, j, k] + Δt *
           ABIII(
                state.G_u[i, j, k, Δt0],
                state.G_u[i, j, k, Δt1],
                state.G_u[i, j, k, Δt2],
        )
    end

    #println("G_u")
    #println(state.G_u[1, 2, 2, Δt0])

    for i=1:env.Nx, j=1:env.Ny+1, k=1:env.Nz_c
        core.v_aux[i, j, k] = state.v_c[i, j, k] + Δt *
            ABIII(
                G_v[i, j, k, Δt0],
                G_v[i, j, k, Δt1],
                G_v[i, j, k, Δt2],
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

    
    Φu     = getSpace!(wksp, :sU)
    Φv     = getSpace!(wksp, :sV)
    DIV_Φu = getSpace!(wksp, :sT)
    DIV_Φv = getSpace!(wksp, :sT)

    calAverage_c2s!(va, core.u_aux, Φu)
    calAverage_c2s!(va, core.v_aux, Φv)
    
    Φu .*= env.Φ_total
    Φv .*= env.Φ_total

    mul2!(DIV_Φu, s_ops.T_DIVx_U,   Φu)
    mul2!(DIV_Φv, s_ops.T_DIVy_V,   Φv)
 
    for i=1:env.Nx, j=1:env.Ny
        core.Φ_aux[i, j] = state.Φ[i, j] - Δt * (
            DIV_Φu[i, j] + DIV_Φv[i, j]
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
 
    for i=1:env.Nx, j=1:env.Ny, k=1:env.Nz_c
        state.u_c[i, j, k] = core.u_aux[i, j, k] - Δt∂Φ∂x[i, j]
        state.v_c[i, j, k] = core.v_aux[i, j, k] - Δt∂Φ∂y[i, j]
    end
 
    #projVertical_c2f!(core.va, state.u_c, state.u_f)
    #projVertical_c2f!(core.va, state.v_c, state.v_f)
   
end

