
include("../../share/constants.jl")

mutable struct DynSlave

    model        :: Dyn.DynModel

    ocn_env      :: OcnEnv
    shared_data  :: SharedData

    function DynSlave(
        ocn_env      :: OcnEnv,
        shared_data  :: SharedData,
    )
       
        mi = ModelMap.MapInfo{Float64}(ocn_env.hrgrid_file)
        gi = PolelikeCoordinate.CurvilinearSphericalGridInfo(;
            R=Re,
            Ω=Ωe,
            Nx=mi.nx,
            Ny=mi.ny,
            c_lon=mi.xc,
            c_lat=mi.yc,
            vs_lon=mi.xv,
            vs_lat=mi.yv,
            area=mi.area,
            angle_unit=:deg,
        )
        
        model = Dyn.DynModel(
            gi      = gi,
            Δt      = ocn_env.Δt,
            z_bnd_f = ocn_env.z_bnd_f,
            height_level_counts = ocn_env.height_level_counts;
            mask    = ocn_env.mask2_deep
        )         

        return new(
            model, 
            ocn_env,
            shared_data,
        )

    end

end
