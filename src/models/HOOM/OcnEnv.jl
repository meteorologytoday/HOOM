mutable struct OcnEnv
   
    hrgrid_file  :: AbstractString
    topo_file  :: AbstractString
    bg_TS_file :: AbstractString

    Δt :: Float64
    substep_dyn :: Int64
    substep_tmd :: Int64

    Nx :: Int64
    Ny :: Int64


    # Fine grids is the total ocean grids
    # tcr and mld cores use fine grids
    # dyn uses coarse grid.  
    Nz_f :: Int64           # number of fine   grids in z direction ( tcr, mld )
    Nz_c :: Int64           # number of coarse grids in z direction ( dyn )
    
    z_bnd_f :: AbstractArray{Float64, 1}           # z boundaries of fine grid
    height_level_counts :: AbstractArray{Int64, 1} # how many layers in fine grid is a coarse grid layer
    
    # First two tracers are T and S. Passive tracers starts from 3.
    NX         :: Int64              # Number of tracers
    NX_passive :: Int64              # Number of passive tracers

    deep_threshold :: Float64                  # threshold where the ocean column is 'deep'

    Kh_m :: Float64   # Horizontal diffusion coe of momentum.
    Kv_m :: Float64   # Vertical   diffusion coe of momentum.

    Kh_X :: AbstractArray{Float64, 1}  # Horizontal diffusion coe of tracers
    Kv_X :: AbstractArray{Float64, 1}  # Vertical   diffusion coe of tracers

    # mld_min, mld_max here are just one numbers.
    # when used in mld core, they will be truncated according to each grid's topography
    mld_min :: Float64    # minimum value of mixed-layer depth
    mld_max :: Float64    # maximum value of mixed-layer depth
    
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



    # Technically mask2 == mask3[:, :, 1]
    mask2      :: Union{AbstractArray{Float64, 2}, Nothing} # where surface ocean grid is active 
    mask2_deep :: Union{AbstractArray{Float64, 2}, Nothing} # where ocean is deep
    topo       :: Union{AbstractArray{Float64, 2}, Nothing} # z value of topography. Negative is below sea surface.
    
    #mask3      :: AbstractArray{Float64, 3}     # where ocean grid is active in 3D grid

    function OcnEnv(
        hrgrid_file :: String,
        topo_file :: String,
        bg_TS_file :: String,
        Δt :: Float64,
        substep_dyn,
        substep_tmd,
        Nz_f :: Int64,
        Nz_c :: Int64,
        z_bnd_f :: AbstractArray{Float64, 1},
        height_level_counts :: AbstractArray{Int64, 1},
        NX_passive :: Int64,
        deep_threshold :: Float64,
        Kh_m  :: Float64,
        Kv_m  :: Float64,
        Kh_X  :: AbstractArray{Float64, 1},
        Kv_X  :: AbstractArray{Float64, 1},
        mld_min :: Float64,
        mld_max :: Float64,
        R       :: Float64,
        ζ       :: Float64,
    )

        local topo

        mi = ModelMap.MapInfo{Float64}(hrgrid_file)
        Dataset(topo_file, "r") do ds
            topo = - (ds["depth"][:] |> nomissing)
        end

        Nx = mi.nx
        Ny = mi.ny
        mask2      = copy(mi.mask)
        
        if size(mask2) != (Nx, Ny)
            throw(ErrorException("Size of mask does not match ocn env."))
        end

        if size(topo) != size(mask2)
            throw(ErrorException("Size of topo and mask does not match"))
        end


        mask2_deep = copy(mask2)
        mask2_deep[topo .> - deep_threshold] .= 0.0

        NX = NX_passive + 2   # 1 = T, 2 = S 

        return new(
            hrgrid_file,
            topo_file,
            bg_TS_file,
            Δt,
            substep_dyn,
            substep_tmd,
            Nx,
            Ny,
            Nz_f,
            Nz_c,
            z_bnd_f,
            height_level_counts,
            NX,
            NX_passive,
            deep_threshold,
            Kh_m,
            Kv_m,
            Kh_X,
            Kv_X,
            mld_min,
            mld_max,
            R,
            ζ,
            mask2,
            mask2_deep,
            topo,
        )

    end

end

function saveOcnEnv(
    env      :: OcnEnv,
    filename :: String,
)
    Dataset(filename, "c") do ds

        defDim(ds, "Nx", env.Nx)
        defDim(ds, "Ny", env.Ny)
        defDim(ds, "Nz_c",  env.Nz_c)
        defDim(ds, "Nz_f",  env.Nz_f)
        defDim(ds, "Nz_fW", env.Nz_f+1)
        defDim(ds, "NX", env.NX)

        ds.attrib["hrgrid_file"] = env.hrgrid_file
        ds.attrib["topo_file"] = env.topo_file
        ds.attrib["bg_TS_file"] = env.bg_TS_file

        ds.attrib["Nz_f"] = env.Nz_f
        ds.attrib["Nz_c"] = env.Nz_c

        ds.attrib["NX"] = env.NX
        ds.attrib["NX_passive"] = env.NX_passive
        ds.attrib["deep_threshold"] = env.deep_threshold
        ds.attrib["Kh_m"] = env.Kh_m
        ds.attrib["Kv_m"] = env.Kv_m
        ds.attrib["mld_min"] = env.mld_min
        ds.attrib["mld_max"] = env.mld_max
        ds.attrib["R"] = env.R
        ds.attrib["zeta"] = env.zeta

        
        for (varname, vardata, vardim, attrib) in [
            ("z_bnd_f",             env.z_bnd_f,             ("Nz_fW",), Dict()),
            ("height_level_counts", env.height_level_counts, ("Nz_c",),  Dict()),
            ("Kh_X",                env.Kh_X,                ("NX",),  Dict()),
            ("Kv_X",                env.Kv_X,                ("NX",),  Dict()),
        ] 

            if ! haskey(ds, varname)
                var = defVar(ds, varname, Float64, vardim)
                var.attrib["_FillValue"] = 1e20
            end

            var = ds[varname]
            
            for (k, v) in attrib
                var.attrib[k] = v
            end

            rng = []
            for i in 1:length(vardim)-1
                push!(rng, Colon())
            end
            push!(rng, 1:size(vardata)[end])
            var[rng...] = vardata

        end
    end
end

function loadOcnEnv(
    filename :: String,
)

    local env

    Dataset(filename, "r") do ds
        env = OcnEnv(
            ds.attrib["hrgrid_file"],
            ds.attrib["topo_file"],
            ds.attrib["bg_TS_file"],
            ds.attrib["Nz_f"],
            ds.attrib["Nz_c"],
            ds["z_bnd_f"][:] |> nomissing,
            ds["height_level_counts"][:] |> nomissing,
            ds.attrib["NX"],
            ds.attrib["NX_passive"],
            ds.attrib["deep_threshold"],
            ds.attrib["Kh_m"],
            ds.attrib["Kv_m"],
            ds["Kh_X"][:] |> nomissing,
            ds["Kv_X"][:] |> nomissing,
            ds.attrib["mld_min"],
            ds.attrib["mld_max"],
            ds.attrib["R"],
            ds.attrib["zeta"],
        )
    end

end

