mutable struct TmdEnv
    
    gi         :: Union{PolelikeCoordinate.GridInfo, Nothing}

    Δt         :: Float64
    Δt_substep :: Float64
    substeps   :: Int64

    Nx         :: Int64           # Number of columns in i direction
    Ny         :: Int64           # Number of columns in j direction
    Nz         :: Int64           # Number of layers  in k direction
    NX         :: Int64   
    NX_passive :: Int64
 
    z_bnd  :: AbstractArray{Float64, 1} # Unmasked z_bnd_av bone
    topo     :: AbstractArray{Float64, 2} # Depth of the topography. Negative value if it is underwater
    z_bnd_av   :: AbstractArray{Float64, 3} # Actuall z_bnd_av coordinate masked by topo
    Nz_av       :: AbstractArray{Int64, 2} # Number of layers that is active

    Kh_X      :: AbstractArray{Float64, 1}           # Horizontal diffusion coe of temperature
    Kv_X      :: AbstractArray{Float64, 1}           # Vertical   diffusion coe of temperature

    mask2    :: AbstractArray{Float64, 2}
    mask3    :: AbstractArray{Float64, 3}
    noflux_x_mask3  :: AbstractArray{Float64, 3}
    noflux_y_mask3  :: AbstractArray{Float64, 3}

    h_ML_min :: AbstractArray{Float64, 2}
    h_ML_max :: AbstractArray{Float64, 2}
    we_max   :: Float64

    # Radiation Scheme
    # The parameterization is referenced to Paulson and Simpson (1977).
    # I made assumption that ζ1→0.0 as adapted in Oberhuber (1993).
    # This means `R` portion of the irradiance is going to be treated
    # as surface heat fluxes (δ like absorption).
    # 
    # Reference: 
    #
    # 1. Oberhuber, J. M. (1993). Simulation of the Atlantic circulation with a coupled sea
    #    ice-mixed layer-isopycnal general circulation model. Part I: Model description.
    #    Journal of Physical Oceanography, 23(5), 808-829.
    #
    # 2. Paulson, C. A., & Simpson, J. J. (1977). Irradiance measurements in the upper ocean.
    #    Journal of Physical Oceanography, 7(6), 952-956.
    #
    R               :: Float64   # Fast absorption portion of sunlight.
    ζ               :: Float64   # Light penetration depth of DO ( = ζ2 in Paulson and Simpson (1977) )

    t_X_wr :: AbstractArray{Float64, 1}
    X_wr   :: AbstractArray{Float64, 4}
    
    MLT_scheme            :: Symbol
    radiation_scheme      :: Symbol
    convective_adjustment :: Bool
    use_Qflux             :: Bool
    finding_Qflux         :: Bool 
 
    function TmdEnv(;
        gi         :: PolelikeCoordinate.GridInfo,
        Δt         :: Float64,
        substeps   :: Int64,
        z_bnd      :: AbstractArray{Float64, 1},
        topo       :: Union{AbstractArray{Float64, 2}, Nothing},
        Kh_X       :: AbstractArray{Float64, 1},
        Kv_X       :: AbstractArray{Float64, 1},
        we_max     :: Float64,
        R          :: Float64,  # See Paulson and Simpson (1977) Type I clear water
        ζ          :: Float64,  # See Paulson and Simpson (1977) Type I clear water
        MLT_rng    :: AbstractArray{Float64, 1},
        t_X_wr     :: Union{AbstractArray{Float64, 1}, Nothing},
        X_wr       :: Union{AbstractArray{Float64, 4}, Nothing},
        NX_passive :: Int64,
        mask2      :: Union{AbstractArray{Float64, 2}, Nothing},
        MLT_scheme    :: Symbol,
        radiation_scheme :: Symbol,
        convective_adjustment :: Bool,
        use_Qflux     :: Bool,
        finding_Qflux :: Bool,
        verbose    :: Bool = false,
    )
       
        Nx = gi.Nx
        Ny = gi.Ny

        if mask2 == nothing
            mask2 = ones(Float64, Nx, Ny)
        end
 
        mask_idx = (mask2 .== 1.0)

        # Arrage like (2, cnt) instead of (cnt, 2) to
        # enhance speed through memory cache
        valid_idx = zeros(Int64, 2, sum(mask_idx))
        
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

        # ===== [BEGIN] topo, mask, h_ML_min, h_ML_max =====
        # Min/max of ML is tricky because it cannot be
        # deeper than the bottom boundary
        # Also, in real data topo can be 0 and not masked out
      
        if topo == nothing
            topo = zeros(Float64, Nx, Ny)
            topo .= z_bnd[end]
        end
        # ===== [END] topo, mask, h_ML_min, h_ML_max =====

        # ===== [BEGIN] z coordinate =====
        z_bnd = copy(z_bnd)
        Nz = length(z_bnd) - 1

        Nz_av    = zeros(Int64,           Nx, Ny)
        z_bnd_av = zeros(Float64, Nz + 1, Nx, Ny)

        z_bnd_av  .= NaN
        Nz_av     .= 0

        for i=1:Nx, j=1:Ny

            if mask2[i, j] == 0
                continue
            end

            # Determine Nz_av

            # Default is that topo is deeper than
            # the bottom of z_bnd
            _Nz_av = Nz
            for k=2:length(z_bnd)
                if z_bnd[k] <= topo[i, j]
                    _Nz_av = k-1
#                    println(format("This topo gets: z_bnd[{:d}] = {:f}, _topo[{:d},{:d}]={:f}", k, z_bnd[k], i, j, topo[i,j]))
                    break
                end
            end

            Nz_av[i, j] = _Nz_av

            # Construct vertical coordinate
            z_bnd_av[1:_Nz_av, i, j] = z_bnd[1:_Nz_av]

            z_bnd_av[_Nz_av+1, i, j] = max(topo[i, j], z_bnd[_Nz_av+1])

        end
        
        # ===== [END] z coordinate =====

        # ===== [BEGIN] construct h_ML limits =====
        h_ML_min    = zeros(Float64, Nx, Ny)
        h_ML_max    = zeros(Float64, Nx, Ny)

        if MLT_rng[2] <= MLT_rng[1]
            throw(ErrorException("MLT_rng should be 2-element with second one larger than the first one."))
        end

        h_ML_min .= MLT_rng[1]
        h_ML_max .= MLT_rng[2]


        # Detect and fix h_ML_{max,min}
        fitMLTToTopo!(
            h_ML_max  = h_ML_max,
            h_ML_min  = h_ML_min,
            topo      = topo,
            mask2     = mask2,
            coord_min = z_bnd[end],
        )
        
        # ===== [END] construct h_ML limits =====

        NX = NX_passive + 2   # T, S, b


        if X_wr == nothing

            X_wr = zeros(Float64, Nz, Nx, Ny, NX)

        elseif size(X_wr) != (Nz, Nx, Ny, NX)
                throw(ErrorException("Invalid size of weak restoration profiles." * string(size(X_wr))))

        end
            
        if t_X_wr == nothing
            t_X_wr = zeros(Float64, NX)
            t_X_wr .= NaN
        end

        # ===== [END] Climatology =====


        # ===== [BEGIN] Mask out data =====

        mask3           = zeros(Float64, Nz, Nx, Ny)
        noflux_x_mask3  = zeros(Float64, Nz, Nx, Ny)
        noflux_y_mask3  = zeros(Float64, Nz, Nx, Ny+1)
        mask3 .= 0.0

        for i=1:Nx, j=1:Ny
            if mask2[i, j] == 1.0
                mask3[1:Nz_av[i, j], i, j] .= 1.0
            end
        end

        # no flow try to avoid the horizontal flow when topography is cutting through a grid box.
        # x
        for i=1:Nx, j=1:Ny
            i_w = cyc(i  , Nx)
            i_e = cyc(i-1, Nx)
            for k=1:Nz_av[i, j]
                noflux_x_mask3[k, i, j] = ( 
                    mask3[k, i_w, j] == 0.0
                 || mask3[k, i_e, j] == 0.0
                 || k >= Nz_av[i_w, j]
                 || k >= Nz_av[i_e, j]
                ) ? 0.0 : 1.0
            end
        end
        
        # y
        for i=1:Nx, j=2:Ny-1     # j=1 and Ny+1 are 0 because these are the boundaries. Need to be careful doing parallization.
            j_s = j
            j_n = j+1
            for k=1:Nz_av[i, j]
                noflux_y_mask3[k, i, j] = (
                    mask3[k, i, j_s] == 0.0
                 || mask3[k, i, j_n] == 0.0
                 || k >= Nz_av[i, j_s]
                 || k >= Nz_av[i, j_n]
                ) ? 0.0 : 1.0
            end
        end


        # ===== [BEGIN] check scheme =====
        if ! ( radiation_scheme in [:step, :exponential_decay] )
                ErrorException("radiation_scheme only accept: `:step` and `:exponential_decay`. I got " * string(radiation_scheme)) |> throw
        end
        # ===== [END] check scheme =====


        # ===== [BEGIN] check integrity =====
        # Check topography, h_ML_min/max and z_bnd_av
        for i=1:Nx, j=1:Ny
            if mask2[i, j] == 0
                continue
            end

            if ! (0 >= - h_ML_min[i, j] >= - h_ML_max[i, j] >= z_bnd_av[Nz_av[i, j] + 1, i, j] >= topo[i, j])
                println("idx: (", i, ", ", j, ")")
                println("h_ML_min: ", h_ML_min[i, j])
                println("h_ML_max: ", h_ML_max[i, j])
                println("z_deepest: ", z_bnd_av[Nz_av[i, j] + 1, i, j])
                println("topo: ", topo[i, j])
                ErrorException("Relative relation is wrong") |> throw
            end
        end

        
        # Check if there is any hole in climatology 
       
        for x=1:NX
            if ! isnan(t_X_wr[x])
                checkDataHoles3(
                    mask3   = mask3,
                    data    = view(X_wr, :, :, :, x),
                    varname = format("{:02d}", x),
                )
            end
        end

        # ===== [END] check integrity =====

        ocn = new(
            gi,
            Δt, Δt/substeps, substeps,
            Nx, Ny, Nz, NX, NX_passive,
            z_bnd, topo,
            z_bnd_av, Nz_av,

            Kh_X, Kv_X,
            mask2, mask3,
            noflux_x_mask3, noflux_y_mask3,
            h_ML_min, h_ML_max, we_max,

            R, ζ,
            
            t_X_wr, X_wr,
           
            MLT_scheme, 
            radiation_scheme,
            convective_adjustment,
            use_Qflux,
            finding_Qflux,
        )
    end
end


function fitMLTToTopo!(;
    h_ML_max  :: AbstractArray{Float64, 2},
    h_ML_min  :: AbstractArray{Float64, 2},
    topo      :: AbstractArray{Float64, 2},
    mask2     :: AbstractArray{Float64, 2},
    coord_min :: Float64,                      # deepest z-coord the model is capable of
    verbose   :: Bool = false
)

    Nx, Ny = size(topo)
    max_thickness = - coord_min

    for i=1:Nx, j=1:Ny

        if mask2[i, j] == 0
            h_ML_max[i, j] = NaN
            h_ML_min[i, j] = NaN
            continue
        end

        hmax = h_ML_max[i, j]
        hmin = h_ML_min[i, j]
        hbot = - topo[i, j]

        if hbot < 0
            throw(ErrorException(format("Topography is above sea-surface and is not masked at idx ({:d}, {:d})", i, j)))
        elseif hbot == 0
            #throw(ErrorException(format("Topography is zero at idx ({:d}, {:d})", i, j)))
            verbose && println(format("Topography is zero at idx ({:d}, {:d})", i, j))
        end

        if hmax < hmin
            throw(ErrorException(format("h_ML_max must ≥ h_ML_min. Problem happens at idx ({:d}, {:d})", i, j)))
        end

        if hmin > max_thickness
            verbose && println(format("Point ({},{}) got h_min {:.2f} which is larger than max_thickness {}. Tune h_ML_min to match it.", i, j, hmin, coord_max))
            hmin = max_thickness
        end

        if hmax > max_thickness
            verbose && println(format("Point ({},{}) got h_max {:.2f} which is larger than max_thickness {}. Tune h_ML_max to match it.", i, j, hmax, coord_max))
            hmax = max_thickness
        end

        if hmin > hbot
            verbose && println(format("Point ({},{}) got depth {:.2f} which is smaller than h_ML_min {}. Tune h_ML_min/max to depth.", i, j, hbot, hmin))
            #hbot = hmin
            hmin = hbot
        end

        if hmax > hbot
            verbose && println(format("Point ({},{}) got depth {:.2f} which is smaller than h_ML_max {}. Tune the h_ML_max to depth.", i, j, hbot, hmax))
            hmax = hbot
        end

        h_ML_min[i, j] = hmin
        h_ML_max[i, j] = hmax
        topo[i, j]     = -hbot

    end
  
end




function checkDataHoles3(;
    data  :: AbstractArray{Float64, 3},
    mask3 :: AbstractArray{Float64, 3},
    varname = "UNKNOWN"
)

    Nz, Nx, Ny = size(mask3)
    valid_ocn_idx  = (mask3 .== 1)
    valid_grids    = sum(valid_ocn_idx)

    total_grid_pts = reduce(*, size(mask3))
    
    finite_data = sum(isfinite.(data[valid_ocn_idx]))

    if finite_data < valid_grids
        throw(ErrorException(format("Data {} has {} holes!", varname, valid_grids - finite_data)))
    end

    
end
