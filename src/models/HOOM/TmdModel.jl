mutable struct TmdModel

    env   :: TmdEnv
    state :: TmdState
    core  :: TmdCore

    function TmdModel(;
        gi,
        Δt,
        z_bnd,
        topo        = nothing,
        Kh_X        = [1e3,  1e3 ],
        Kv_X        = [1e-5, 1e-5],
        we_max      = 1e-2,
        R           = 0.58,  # See Paulson and Simpson (1977) Type I clear water
        ζ           = 23.0,  # See Paulson and Simpson (1977) Type I clear water
        MLT_rng     = [10.0, 1000.0],
        t_wr_X      = nothing,
        X_wr        = nothing,
        NX_passive  = 0,
        mask        = nothing,
    )

        env = TmdEnv(;
            gi         = gi,
            Δt         = Δt,
            z_bnd      = z_bnd,
            topo       = topo,
            Kh_X       = Kh_X,
            Kv_X       = Kv_X,
            we_max     = we_max,
            R          = R,
            ζ          = ζ,
            MLT_rng    = MLT_rng, 
            NX_passive = NX_passive,
            mask       = mask,
        )

        state = TmdState(env)
        core  = TmdCore(env, state)

        return(
            env, state, core
        )

    end
end
