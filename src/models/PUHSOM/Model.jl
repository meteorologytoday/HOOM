mutable struct Model

    env            :: OcnEnv
    shared_data    :: SharedData
    job_dist_info  :: JobDistributionInfo
    data_exchanger :: DataExchanger

#    function Model()
#        return new()
#    end
end

