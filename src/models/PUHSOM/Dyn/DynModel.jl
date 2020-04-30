mutable struct DynModel

    env     :: DynEnv
    state   :: DynState
    core    :: DynCore
    
    
    function DynModel(;
        gi :: PolelikeCoordinate.GridInfo,
        Δt,
        Dh,
        z_bnd_f, height_level_counts,
        NX_passive=0,
        mask=nothing,
    )

        env = DynEnv(;
            gi = gi,
            Δt = Δt,
            Dh = Dh,
            Nx = gi.Nx,
            Ny = gi.Ny,
            z_bnd_f = z_bnd_f, 
            height_level_counts = height_level_counts,
            NX_passive = NX_passive,
            mask=mask,
        )

        state = DynState(env, :local)
        core  = DynCore(env, state)

        #=
        state.u_c .= 1000
        state.v_c .= 1000

        uu = copy(state.u_c)
        vv = copy(state.v_c)

        mul3!(uu, core.s_ops.filter_U, state.u_c)
        mul3!(vv, core.s_ops.filter_V, state.v_c)
        
        state.u_c .-= uu
        state.v_c .-= vv
       
        state.u_c[uu .== 0] .= NaN 
        state.v_c[vv .== 0] .= NaN 
        =#


        return new(
            env,
            state,
            core,
        )

    end
end


