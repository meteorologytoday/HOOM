mutable struct Cols
    z_bnd_av      :: Array{SubArray}
    T             :: Array{SubArray}
    S             :: Array{SubArray}
    b             :: Array{SubArray}
    dz_W          :: Array{SubArray}
    rad_decay_coe :: Array{SubArray}
    rad_absorp_coe:: Array{SubArray}
    function Cols()
        return new()
    end
end



mutable struct TmdCore    # Adam Bashford

    cols           :: Cols

    dz_W           :: AbstractArray{Float64, 3}   
    dz_T           :: AbstractArray{Float64, 3}   

    rad_decay_coe  :: AbstractArray{Float64, 3}
    rad_absorp_coe :: AbstractArray{Float64, 3}

    ASUM           :: Any#AdvectionSpeedUpMatrix

    wksp           :: Workspace

    XFLUX_top    :: AbstractArray{Float64, 3}
    XFLUX_bot    :: AbstractArray{Float64, 3}
    XFLUX_CONV   :: AbstractArray{Float64, 4}
    XFLUX_CONV_h :: AbstractArray{Float64, 4}
    XFLUX_DEN_x  :: AbstractArray{Float64, 4}
    XFLUX_DEN_y  :: AbstractArray{Float64, 4}
    XFLUX_DEN_z  :: AbstractArray{Float64, 4}

    acc_vars     :: AccumulativeVariables
    
    valid_idx    :: AbstractArray{Int64, 2}

    current_substep :: Int64

    function TmdCore(env, state, diag, forcing)
        
        Nx = env.Nx
        Ny = env.Ny
        Nz = env.Nz
        NX = env.NX

        mask2 = env.mask2

        dz_W     = zeros(Float64, Nz    , Nx, Ny)
        dz_T     = zeros(Float64, Nz - 1, Nx, Ny)

        dz_W      .= NaN
        dz_T      .= NaN
        for i=1:Nx, j=1:Ny
            if mask2[i, j] == 0
                continue
            end
            _Nz_av = env.Nz_av[i, j]
            dz_W[ 1:_Nz_av,   i, j] = env.z_bnd_av[1:_Nz_av, i, j] - env.z_bnd_av[2:_Nz_av+1, i, j]
            dz_T[1:_Nz_av-1, i, j] = (dz_W[1:_Nz_av-1, i, j] + dz_W[2:_Nz_av, i, j]) / 2.0
        end
        
        rad_decay_coe, rad_absorp_coe = genRadCoe(
            Nx       = Nx,
            Ny       = Ny,
            Nz       = Nz,
            Nz_av    = env.Nz_av,
            mask2    = mask2,
            ζ        = env.ζ,
            dz_W     = dz_W,
            z_bnd_av = env.z_bnd_av,
        )

        # making columnwise views
        cols = Cols()
        for (var, ref) in Dict(
            :z_bnd_av        => env.z_bnd_av,
            :T               => state.T,
            :S               => state.S,
            :b               => state.b,
            :dz_W            => dz_W,
            :rad_decay_coe   => rad_decay_coe,
            :rad_absorp_coe  => rad_absorp_coe,
        )
            setfield!(cols, var, genColView(ref))
        end

        println("Making ASUM")
        @time ASUM = nothing; 
        #=AdvectionSpeedUpMatrix(;
                gi = env.gi,
                Nx = Nx,
                Ny = Ny,
                Nz_bone = Nz,
                Nz = env.Nz_av,
                mask3 = env.mask3,
                noflux_x_mask3 = env.noflux_x_mask3,
                noflux_y_mask3 = env.noflux_y_mask3,
                Δzs = dz_T,
                hs  = dz_W,
        )=#

        wksp = Workspace(;
            Nx = Nx,
            Ny = Ny,
            Nz = Nz,
            T = 8,
            U = 1,
            V = 1,
            W = 1,
            sT = 0,
            sU = 0,
            sV = 0,
        ) 


        XFLUX_top    = zeros(Float64, Nx, Ny, NX)
        XFLUX_bot    = zeros(Float64, Nx, Ny, NX)
        XFLUX_CONV   = zeros(Float64, Nz, Nx, Ny, NX)
        XFLUX_CONV_h = zeros(Float64, Nz, Nx, Ny, NX)
        XFLUX_DEN_x  = zeros(Float64, Nz, Nx, Ny, NX)
        XFLUX_DEN_y  = zeros(Float64, Nz, Nx, Ny+1, NX)
        XFLUX_DEN_z  = zeros(Float64, Nz+1, Nx, Ny, NX)
   
        acc_vars = AccumulativeVariables(Nx, Ny, Nz, NX)
 

        valid_idx = zeros(Int64, 2, sum(mask2 .== 1.0))
        let k = 1
            for idx in CartesianIndices((Nx, Ny))
                if mask2[idx] == 1.0
                    valid_idx[1, k] = idx[1]
                    valid_idx[2, k] = idx[2]

                    k += 1
                end
            end

            if k != size(valid_idx)[2] + 1
                throw(ErrorException("Initialization error making `valid_idx`"))
            end
        end

        current_substep = 1
    
        new(


            cols,
            dz_W,
            dz_T,
            rad_decay_coe,
            rad_absorp_coe,
            ASUM,
            wksp,         
            XFLUX_top,
            XFLUX_bot,
            XFLUX_CONV,
            XFLUX_CONV_h,
            XFLUX_DEN_x,
            XFLUX_DEN_y,
            XFLUX_DEN_z,
            acc_vars,
            valid_idx,
            current_substep,
        )
    end
end

function genColView(
    arr :: AbstractArray{T, 3}
) where T

    _, Nx, Ny = size(arr)
    view_arr  = Array{SubArray}(undef, Nx, Ny)
    for i=1:Nx, j=1:Ny
        view_arr[i, j] = view(arr, :, i, j)
    end

    return view_arr
end

function genRadCoe(;
    Nx       :: Int64,
    Ny       :: Int64,
    Nz       :: Int64,
    Nz_av    :: AbstractArray{Int64,   2},
    mask2    :: AbstractArray{Float64, 2},
    ζ        :: Float64,
    dz_W     :: AbstractArray{Float64, 3},
    z_bnd_av :: AbstractArray{Float64, 3},
)

    rad_decay_coe  = zeros(Float64, Nz, Nx, Ny)
    rad_absorp_coe = zeros(Float64, Nz, Nx, Ny)

    for i=1:Nx, j=1:Ny
        
        if mask2[i, j] == 0.0
            continue
        end

        for k=1:Nz_av[i, j]
            rad_decay_coe[k, i, j]  = exp(z_bnd_av[k, i, j] / ζ)         # From surface to top of the layer
            rad_absorp_coe[k, i, j] = 1.0 - exp(- dz_W[k, i, j] / ζ)
        end
        
        # Since we assume the bottome of ocean absorbs anything
        rad_absorp_coe[Nz_av[i, j], i, j] = 1.0
    end

    return rad_decay_coe, rad_absorp_coe
        
end
