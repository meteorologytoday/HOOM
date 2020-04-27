
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

        data_exchanger = DataExchanger([
            :MOM,
            :THERMO,
        ])

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
    du_there = sd.data_units

    bindings = (
        ([:MOM],    DataUnit(:u_c, :cU, :xyz, s.u_c, false), :u_c),
        ([:MOM],    DataUnit(:v_c, :cV, :xyz, s.v_c, false), :v_c),
        ([:THERMO], DataUnit(:b_c, :cT, :xyz, s.b_c, false), :b_c),
        ([:THERMO], DataUnit(:Φ  , :sT, :xy , s.Φ,   false), :Φ  ),
    )

    for (group_labels, here, there_key) in bindings
       
        println("Doing : ", here.id, "; ", du_there[there_key].id) 

        createBinding!(
            de,
            here,
            du_there[there_key],
            Colon(),
            Colon(),
            labels = group_labels,
        )
    end
end
