mutable struct OcnEnv
   
    grid_file  :: AbstractString
    topo_file  :: AbstractString
    bg_TS_file :: AbstractString

    Δt :: Float64
    substep_dyn :: Int64
    substep_tcr :: Int64
    substep_mld :: Int64

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
    NX   :: Int64                         # Number of tracers
    NX_passive :: Int64                   # Number of passive tracers

    deep_threshold :: Float64                  # threshold where the ocean column is 'deep'

    Kh_m :: Float64   # Horizontal diffusion coe of momentum.
    Kv_m :: Float64   # Vertical   diffusion coe of momentum.

    Kh_X :: AbstractArray{Float64, 2}  # Horizontal diffusion coe of tracers
    Kv_X :: AbstractArray{Float64, 2}  # Vertical   diffusion coe of tracers

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



    #= these variable should be instatiated in cores
    # Technically mask2 == mask3[:, :, 1]
    mask2      :: AbstractArray{Float64, 2}     # where surface ocean grid is active 
    mask3      :: AbstractArray{Float64, 3}     # where ocean grid is active in 3D grid

    mask2_deep :: AbstractArray{Float64, 2}     # where ocean is deep
    topo       :: AbstractArray{Float64, 2}     # z value of topography. Negative is below sea surface.
    =#

end


