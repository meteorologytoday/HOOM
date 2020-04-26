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
        NX_passive  = 0,
        t_X_wr      = nothing,
        X_wr        = nothing,
        mask2       = nothing,
        radiation_scheme = :exponential_decay,
        convective_adjustment = true,
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
            t_X_wr     = t_X_wr,
            X_wr       = X_wr,
            mask2      = mask2,
            radiation_scheme = radiation_scheme,
            convective_adjustment = convective_adjustment,
        )

        state = TmdState(env)
        core  = TmdCore(env, state)

        return new(
            env, state, core
        )

    end
end