#include("../../share/constants.jl")

mutable struct DynSlave

    model        :: Dyn.DynModel

    ocn_env      :: OcnEnv
    shared_data  :: SharedData

    data_exchanger :: DataExchanger

    function DynSlave(
        ocn_env      :: OcnEnv,
        shared_data  :: SharedData,
    )
       
        gi = PolelikeCoordinate.genGridInfo(
            ocn_env.hrgrid,
        )
        
        model = Dyn.DynModel(
            gi      = gi,
            Δt      = ocn_env.Δt / ocn_env.substeps_dyn,
            Dh      = ocn_env.Kh_m,
            z_bnd_f = ocn_env.z_bnd_f,
            height_level_counts = ocn_env.height_level_counts;
            mask    = ocn_env.mask2_deep
        )         

        data_exchanger = DataExchanger([
            :FR_TMD,
            :TO_TMD,
            :TO_MAS,
            :TEST,
        ])

        return new(
            model, 
            ocn_env,
            shared_data,
            data_exchanger,
        )

    end

end

