mutable struct TmdModel

    env     :: TmdEnv
    state   :: TmdState
    core    :: TmdCore
    
    
    function TmdModel(;
        gi,
        Δt,
        z_bnd,
        NX_passive=0,
        mask=nothing,
    )

        env = TmdEnv(;
            Δt         = Δt,
            gi         = gi,
            z_bnd      = z_bnd, 
            NX_passive = NX_passive,
            mask       = nothing,
        )

        state = DynState(env, :local)
        core  = DynCore(env, state)

        return new(
            env,
            state,
            core,
        )

    end
end


