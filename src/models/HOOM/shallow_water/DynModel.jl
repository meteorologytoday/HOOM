mutable struct DynModel

    env     :: DynEnv
    state   :: DynState
    core    :: DynCore
    
    
    function DynModel(;
        gi,
        Δt,
        z_bnd_f, height_level_counts,
        NX_passive=0,
        mask=nothing,
    )

      
        env = DynEnv(;
            Δt = Δt,
            gi = gi,
            Nx = gi.Nx,
            Ny = gi.Ny,
            z_bnd_f = z_bnd_f, 
            height_level_counts = height_level_counts,
            NX_passive = NX_passive,
            mask=nothing,
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


