function getCoordinateVariable(m::Model)
    return Dict(
#        "f"               => ( m.env.gi.c_f,                                                ("Nx", "Ny",) ),
    )
end

function getCompleteVariableList(m::Model)
        d = m.shared_data.data_units
        return Dict(
            "Phi"               => ( d[:Φ].data,                               ("Nx", "Ny",) ),
        )
end


function getBasicRecorder(
    m :: Model,
)
    

    # Setting up recorder
    complete_varlist = getCompleteVariableList(m)
    varlist = []
    for varname in keys(complete_varlist)
        println(format("Using varaible: {:s}", varname))
        push!(varlist, ( varname, complete_varlist[varname]... ) )
    end

    coord_varlist = getCoordinateVariable(m)
    cvarlist = []
    for varname in keys(coord_varlist)
        println(format("Using varaible: {:s}", varname))
        push!(cvarlist, ( varname, coord_varlist[varname]... ) )
    end

    return RecordTool.Recorder(
        Dict(
            "NX"      => m.env.NX,
            "Nyp1"    => m.env.Ny+1,
            "Nz_fp1"  => m.env.Nz_f+1,
            "Nx"      => m.env.Nx,
            "Ny"      => m.env.Ny,
            "Nz_f"    => m.env.Nz_f,
            "Nz_c"    => m.env.Nz_c,
            "z_bnd_f" => length(m.env.z_bnd_f),
        ), varlist, Dict(),
        other_vars = cvarlist
    )
     

end

