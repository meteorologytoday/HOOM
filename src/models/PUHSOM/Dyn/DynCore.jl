mutable struct DynCore    # Adam Bashford

    c_ops    :: Union{DynamicAdvSpeedUpMatrix, Nothing}
    s_ops    :: Union{DynamicAdvSpeedUpMatrix, Nothing}
    va       :: Union{VerticalAverager, Nothing}
    Φ_solver :: PhiSolver

    G_idx    :: Dict
    wksp     :: Workspace

    u_aux    :: AbstractArray{Float64, 3}
    v_aux    :: AbstractArray{Float64, 3}
    Φ_aux    :: AbstractArray{Float64, 2}
    
    function DynCore(env, state)
        
        Nx = env.Nx
        Ny = env.Ny
        Nz = env.Nz_c

        va = VerticalAverager(
            z_bnd_f = env.z_bnd_f,
            height_level_counts = env.height_level_counts
        )
        
        println("Making Spatial Operators")
        @time s_ops = DynamicAdvSpeedUpMatrix(;
                gi = env.gi,
                Nz = 1,
                mask2 = env.mask,
        )

        c_ops = s_ops

        #=
        @time c_ops = DynamicAdvSpeedUpMatrix(;
                gi = env.gi,
                Nz = env.Nz_c,
                mask2 = env.mask,
        )
        =#

        Φ_solver = PhiSolver(
            gi    = env.gi,
            mask2 = env.mask,
            α     = 1.0 / (env.Δt^2 * env.Φ_total),
            M     = s_ops,
        )

        #         now   Δt-ago  2Δt-ago
        G_idx = Dict(
            :now           =>  1,
            :one_Δt_ago    =>  2,
            :two_Δt_ago    =>  3,
        )

        println("Creating Workspaces")
        wksp = Workspace(;
            Nx = Nx,
            Ny = Ny,
            Nz_c = env.Nz_c,
            Nz_f = env.Nz_f,
            fT = 5,
            fU = 5,
            fV = 5,
            cT = 5,
            cU = 6,
            cV = 6,
            sT = 5,
            sU = 5,
            sV = 5,
        ) 

        u_aux = zeros(Float64, Nx, Ny,   Nz)
        v_aux = zeros(Float64, Nx, Ny+1, Nz)
        Φ_aux = zeros(Float64, Nx, Ny)

        new(
            c_ops,
            s_ops,
            va,
            Φ_solver,
            G_idx,
            wksp,
            u_aux,
            v_aux,
            Φ_aux,
        )
    end
end


