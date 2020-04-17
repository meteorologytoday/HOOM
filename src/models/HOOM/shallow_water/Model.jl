mutable struct Model

    env     :: Env
    state   :: State
    tcr_adv :: TracerAdv
    dyn_adv :: Any#DynamicAdv

    function Model(;
        gi,
        z_bnd_f, height_level_counts,
        f, ϵ,
        Dh,
        Dv,
        NX_passive = 0,
        Nz_f_av=nothing, H_f=nothing, Δz_f=nothing,
        mask3_f=nothing, noflux_x_mask3_f=nothing, noflux_y_mask3_f=nothing,
    )

      
        datakind=:local 

        env = Env(;
            gi = gi,
            Nx = gi.Nx,
            Ny = gi.Ny,
            z_bnd_f = z_bnd_f, 
            height_level_counts = height_level_counts,
            NX_passive = NX_passive,
            Dh=Dh,
            Dv=Dv,
            Nz_f_av = Nz_f_av,
            H_f = H_f,
            Δz_f = Δz_f,
            mask3_f = mask3_f,
            noflux_x_mask3_f = noflux_x_mask3_f,
            noflux_y_mask3_f = noflux_y_mask3_f,
            datakind=datakind,
            f=f,
            ϵ=ϵ,
        )

        state = State(env, datakind)
        tcr_adv = TracerAdv(env, state, datakind)
        #dyn_adv = DynamicAdv(env, state, datakind)

        return new(
            env,
            state,
            tcr_adv,
            nothing, #dyn_adv,
        )

    end
end


