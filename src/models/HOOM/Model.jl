mutable struct Model
    env
    shared_data
    job_dist_info

    function Model()
        return new()
    end
end

