
include("../../share/constants.jl")

mutable struct DynSlave

    model        :: Dyn.DynModel

    ocn_env      :: OcnEnv
    shared_data  :: SharedData

    data_exchanger :: DataExchanger

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

        data_exchanger = DataExchanger()

        return new(
            model, 
            ocn_env,
            shared_data,
            data_exchanger,
        )

    end

end


function setupBinding!(
    slave :: DynSlave,
)

    de = slave.data_exchanger
    sd = slave.shared_data
    m  = slave.model
    s  = m.state
    du = sd.data_units

    bindings = (
        ("u_c", s.u_c, :xyz, :u_c),
        ("v_c", s.v_c, :xyz, :v_c),
        ("b_c", s.b_c, :xyz, :b_c),
        ("Phi", s.Φ,   :xy , :Phi),
    )

    for (name, data_here, data_here_shape, data_there_key) in bindings
        createBinding!(
            de,
            name, 
            data_here, data_here_shape,
            du[data_there_key].data, du[data_there_key].shape
        )
    end


end
