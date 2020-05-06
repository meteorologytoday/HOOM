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

#=
  0.041038 seconds (106 allocations: 9.164 MiB, 2.54% gc time)
  0.045352 seconds (2.16 M allocations: 32.969 MiB, 2.54% gc time)
  0.026442 seconds (1.44 M allocations: 21.995 MiB, 5.05% gc time)
  0.005663 seconds (21 allocations: 1.572 MiB)
  0.001000 seconds (100 allocations: 3.063 KiB)
=#


    decomposeModes!(model)

#    println("Do Diffusion!")
#    doHDiffusionBarotropic!(model)

#    println("calAuxV!")
    calAuxV!(model)

#    println("calAuxΦ!")
    calAuxΦ!(model)

#    println("solveΦ!")
    solveΦ!(model)

#    model.state.Φ .= model.env.Φ_total

#    println("updateV!")
    updateV!(model)

    #println(format("(u, v) = ({:.2f}, {:.2f})", state.u_total[1, 5, 5], state.v_total[1,5,5]))
end

function doHDiffusionBarotropic!(
    model :: DynModel,
)
 
    state = model.state
    core = model.core
    env = model.env
    wksp = core.wksp
    layers = core.layers 
    #wksp_cU = getSpace!(wksp, :cU)
    #wksp_cV = getSpace!(wksp, :cV)
    #wksp_cU .= state.u_total 
    #wksp_cV .= state.v_total 
    #wksp_cU .= state.u_total 
    #wksp_cV .= state.v_total 
#= 
    for k = 1:env.Nz_c

        solveDiffusion!(
            core.diffusion_solver, :U,
            view(wksp_cU  , :, :, k),
            view(state.u_total, :, :, k),
        )

        solveDiffusion!(
            core.diffusion_solver, :V,
            view(wksp_cV  , :, :, k),
            view(state.v_total, :, :, k),
        )
    end
=#
    wksp_sU = getSpace!(wksp, :sU)
    wksp_sV = getSpace!(wksp, :sV)

    wksp_sU .= state.U
    wksp_sV .= state.V 
 

    solveDiffusion!(
        core.diffusion_solver, :U,
        wksp_sU,
        state.U,
    )

    solveDiffusion!(
        core.diffusion_solver, :V,
        wksp_sV,
        state.V,
    )


    for k = 1:env.Nz_c
        @. layers.u_total[k] = state.U + layers.u[k]
    end
end

function decomposeModes!(
    model :: DynModel
)
 
    state = model.state
    core = model.core
    c_ops = core.c_ops
    s_ops = core.s_ops
    va = core.va
    env = model.env
    wksp = core.wksp
  
    
    # cal barotropic and baroclinic components
    decompose!(total=state.u_total, mean=state.U, anomaly=state.u, va=va)
    decompose!(total=state.v_total, mean=state.V, anomaly=state.v, va=va)

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
  
  
    # cal τx_acc, τy_acc
    τx_acc = getSpace!(wksp, :sU)
    τy_acc = getSpace!(wksp, :sV)

    mul2!(τx_acc, s_ops.U_interp_T, state.τx)
    mul2!(τy_acc, s_ops.V_interp_T, state.τy)
    τx_acc ./= (ρ_fw * env.H_c[1])
    τy_acc ./= (ρ_fw * env.H_c[1])
 
    # cal ∇b
    ∂B∂x   = core.∂B∂x
    ∂B∂y   = core.∂B∂y
    mul3!(core.∂B∂x, c_ops.U_∂x_T, state.B)
    mul3!(core.∂B∂y, c_ops.V_∂y_T, state.B)

    # cal Coriolis force
    fu   = getSpace!(wksp, :cV)
    fv   = getSpace!(wksp, :cU)
    mul3!(fu, c_ops.V_f_U, state.u_total)   # fu on V grid
    mul3!(fv, c_ops.U_f_V, state.v_total)   # fv on U grid
    #println("fu: ", fu[30, 29, 1])   

 
    # ===== [ BEGIN cal (v⋅∇)v ] =====
    #=
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
    mul3!(v∂u∂y, c_ops.U_interp_V, state.v_total)  # store interpolated v into v∂u∂y
    for i=1:env.Nx, j=1:env.Ny, k=1:env.Nz_c
        u∂u∂x[i, j, k]  = state.u[i, j, k] * ∂u∂x[i, j, k]
        v∂u∂y[i, j, k] *=                    ∂u∂y[i, j, k]
    end

    # On V grid
    u∂v∂x = getSpace!(wksp, :cV)
    v∂v∂y = getSpace!(wksp, :cV)
    mul3!(u∂v∂x, c_ops.V_interp_U, state.u_total)  # store interpolated u into u∂v∂x
    for i=1:env.Nx, j=1:env.Ny+1, k=1:env.Nz_c
        u∂v∂x[i, j, k] *=                     ∂v∂x[i, j, k]
        v∂v∂y[i, j, k]  = state.v[i, j, k]  * ∂v∂y[i, j, k]
    end
    =#
    # ===== [ END cal (v⋅∇)v ] =====

    # cal G
    G_idx = core.G_idx
    Δt0 = G_idx[:now]
    Δt1 = G_idx[:one_Δt_ago]
    Δt2 = G_idx[:two_Δt_ago]


    G_u = view(state.G_u, :, :, :, Δt0)
    G_v = view(state.G_v, :, :, :, Δt0)

    @. G_u = ∂B∂x + fv
#    G_u .-= u∂u∂x
#    G_u .-= v∂u∂y
 
    @. G_v = ∂B∂y - fu
#    G_v .-= u∂v∂x
#    G_v .-= v∂v∂y


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
#    G_u = state.G_u
#    G_v = state.G_v

    G_u_Δt0 = view(state.G_u, :, :, :, Δt0)
    G_u_Δt1 = view(state.G_u, :, :, :, Δt1)
    G_u_Δt2 = view(state.G_u, :, :, :, Δt2)

    G_v_Δt0 = view(state.G_v, :, :, :, Δt0)
    G_v_Δt1 = view(state.G_v, :, :, :, Δt1)
    G_v_Δt2 = view(state.G_v, :, :, :, Δt2)


    @. core.u_aux = state.u_total + Δt * ABIII(G_u_Δt0, G_u_Δt1, G_u_Δt2)
    @. core.v_aux = state.v_total + Δt * ABIII(G_v_Δt0, G_v_Δt1, G_v_Δt2)

    #=
    @time for i=1:env.Nx, j=1:env.Ny, k=1:env.Nz_c
        core.u_aux[i, j, k] = state.u_total[i, j, k] + Δt *
           ABIII(
                state.G_u[i, j, k, Δt0],
                state.G_u[i, j, k, Δt1],
                state.G_u[i, j, k, Δt2],
           )
    end

    #println("G_u")
    #println(state.G_u[1, 2, 2, Δt0])

    for i=1:env.Nx, j=1:env.Ny+1, k=1:env.Nz_c
        core.v_aux[i, j, k] = state.v_total[i, j, k] + Δt *
            ABIII(
                G_v[i, j, k, Δt0],
                G_v[i, j, k, Δt1],
                G_v[i, j, k, Δt2],
            )
    end
    =#

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

 
    @. core.Φ_aux = state.Φ - Δt * (DIV_Φu + DIV_Φv)
   
    #= 
    for i=1:env.Nx, j=1:env.Ny
        core.Φ_aux[i, j] = state.Φ[i, j] - Δt * (
            DIV_Φu[i, j] + DIV_Φv[i, j]
        )
    end
    =#
end

function solveΦ!(
    model :: DynModel,
)
    env    = model.env
    core   = model.core
    state  = model.state
    solver = core.Φ_solver

    lhs = getSpace!(core.wksp, :sT)
    rhs = getSpace!(core.wksp, :sT)

    lhs_TT = view(getSpace!(core.wksp, :sT), 1:solver.TT_length)
    rhs_TT = view(getSpace!(core.wksp, :sT), 1:solver.TT_length)

    α = solver.α
    
    @. rhs = - core.Φ_aux * α

    mul!(rhs_TT, solver.TT_send_T, view(rhs, :))
 
    ldiv!(
        lhs_TT,
        solver.tool_mtx.MoLap,
        rhs_TT,
    )
        
    mul!(view(state.Φ, :), solver.T_send_TT, lhs_TT)

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

    # TODO: NEED to use AUXILIRAY U, V instead of old U, V 
    ΔtfricU = getSpace!(wksp, :sU)
    ΔtfricV = getSpace!(wksp, :sV)

    ΔtfricU .= state.U
    ΔtfricV .= state.V

    τ = 100*86400.0

    ΔtfricU .*= Δt/τ
    ΔtfricV .*= Δt/τ

    
    for k=1:env.Nz_c
        u_total = view(state.u_total, :, :, k)
        v_total = view(state.v_total, :, :, k)
        u_aux = view(core.u_aux, :, :, k)
        v_aux = view(core.v_aux, :, :, k)

        @. u_total = u_aux - Δt∂Φ∂x - ΔtfricU
        @. v_total = v_aux - Δt∂Φ∂y - ΔtfricV
    end

#= 
    for i=1:env.Nx, j=1:env.Ny, k=1:env.Nz_c
        state.u_total[i, j, k] = core.u_aux[i, j, k] - Δt∂Φ∂x[i, j] - ΔtfricU[i, j]
        state.v_total[i, j, k] = core.v_aux[i, j, k] - Δt∂Φ∂y[i, j] - ΔtfricV[i, j]
    end
 =#
    #projVertical_c2f!(core.va, state.u_total, state.u_f)
    #projVertical_c2f!(core.va, state.v_total, state.v_f)
   
end

