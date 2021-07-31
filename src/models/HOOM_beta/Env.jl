mutable struct Env

    id       :: Integer  # 0 = master, 1, ..., N = workers
    
    gf_filename :: Union{AbstractString, Nothing}
    gf          :: Union{PolelikeCoordinate.GridFile, Nothing}

    gd       :: PolelikeCoordinate.Grid
    gd_slab  :: PolelikeCoordinate.Grid

    z_w      :: AbstractArray{Float64, 1} # Unmasked zs bone

    Dh_T      :: Float64           # Horizontal diffusion coe of temperature
    Dv_T      :: Float64           # Vertical   diffusion coe of temperature
    Dh_S      :: Float64           # Horizontal diffusion coe of salinity
    Dv_S      :: Float64           # Vertical   diffusion coe of salinity

    mask_sT   :: AbstractArray{Float64, 3}
    mask_T    :: AbstractArray{Float64, 3}

    R         :: Float64   # Fast absorption portion of sunlight.
    ζ1        :: Float64   # Light penetration depth of DO ( = ζ2 in Paulson and Simpson (1977) )
    ζ2        :: Float64   # Light penetration depth of DO ( = ζ2 in Paulson and Simpson (1977) )

    τ_TEMP    :: Union{Float64, Nothing}
    τ_SALT    :: Union{Float64, Nothing}

    function Env(;
        id              :: Integer = 0,  
        gf_filename :: AbstractString,
        sub_yrng :: Union{UnitRange, Nothing} = nothing,
        z_w      :: AbstractArray{Float64, 1},
        Dh_T     :: Float64 = 1e3,
        Dv_T     :: Float64 = 1e-5,
        Dh_S     :: Float64 = 1e3,
        Dv_S     :: Float64 = 1e-5,
        τ_TEMP   :: Union{Float64, Nothing},
        τ_SALT   :: Union{Float64, Nothing},
        ϵ        :: Union{AbstractArray{Float64, 2}, Float64, Nothing} = 1.0 / 86400.0,
        qflx_finding_forcing_file :: String = "",
        R        :: Float64 = 0.58,
        ζ1       :: Float64 = 23.0,
        ζ2       :: Float64 = 0.35,
        verbose  :: Bool = false,
    )
        
        # ===== [BEG] GridInfo =====

        # mask =>   lnd = 0, ocn = 1
        gf = PolelikeCoordinate.CurvilinearSphericalGridFile(
            gf_filename;
            R  = 6371229.0,
            Ω  = 2π / (86400 / (1 + 365/365)),
        )

        if id == 0

            gd = nothing
            gd_slab = nothing

        else

            if sub_yrng == nothing
                thorw(ErrorException("Init worker ocean, sub_yrng must be provided."))
            end

            gd      = PolelikeCoordinate.genGrid(gf, z_w ; sub_yrng=sub_yrng) 
            gd_slab = PolelikeCoordinate.genGrid(gf, [0, -1.0]; sub_yrng=sub_yrng) 

        end

        mask_sT = reshape(gf.mask, 1, size(gf.mask)...)
        mask_T  = repeat( mask_sT, outer=(gd.Nz, gd.Nx, gd.Ny) )
        

        # ===== [BEGIN] Q-flux finding forcing file =====
        qflx_finding_cdm = nothing
        
        if id == 0
            # do nothing 
        else
            if sub_yrng == nothing
                thorw(ErrorException("Init worker ocean, sub_yrng must be provided."))
            end

            if qflx_finding_forcing_file != ""

                println("Load qflx file: ", qflx_finding_forcing_file)
                qflx_finding_cdm = CyclicDataManager(;
                    filename     = qflx_finding_forcing_file,
                    varname_time = "time", 
                    varnames     = ["TEMP", "SALT"],
                    beg_time     = 0.0,
                    cyc_time     = 365.0,
                    spatial_rng  = (:, sub_yrng, :),  # in the shape of the file
                    xyz2zxy      = true,
                )
            end
               
        end

        return new(
            id,
            
            gf_filename,
            gf,

            gd,
            gd_slab,
            
            z_w,

            Dh_T,
            Dv_T,
            Dh_S,
            Dv_S,

            mask_sT,
            mask_T,

            R,
            ζ1,
            ζ2,

            τ_TEMP,
            τ_SALT,
        )
    end

end


