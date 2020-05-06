function getCoordinateVariable(m::Model)
    return Dict(
#        "f"               => ( m.env.gi.c_f,                                                ("Nx", "Ny",) ),
    )
end

function getCompleteVariableList(m::Model)
        d = m.shared_data.data_units
        return Dict(
            "T"               => ( view(d[:X].odata, :, :, :, 1),             ("Nx", "Ny", "Nz_f") ),
#            "S"               => ( view(d[:X].odata, :, :, :, 2),             ("Nx", "Ny", "Nz_f") ),
            "swflx"           => ( d[:SWFLX].odata,                           ("Nx", "Ny") ),
            "nswflx"          => ( d[:NSWFLX].odata,                          ("Nx", "Ny") ),
            "X_ML"            => ( view(d[:X_ML].odata, :, :, :, 1),          ("Nx", "Ny", "NX") ),
            "h_ML"            => ( d[:h_ML].odata,                            ("Nx", "Ny") ),
            "Phi"             => ( d[:Φ].odata,                               ("Nx", "Ny",) ),
            "u_total"         => ( d[:u_total_c].odata,                       ("Nx", "Ny",   "Nz_c") ),
            "v_total"         => ( d[:v_total_c].odata,                       ("Nx", "Nyp1", "Nz_c") ),
            "b_ML"            => ( d[:b_ML].odata,                            ("Nx", "Ny") ),
            "b"               => ( d[:b].odata,                               ("Nx", "Ny", "Nz_f") ),
#            "B"               => ( d[:B].odata,                               ("Nx", "Ny", "Nz_f") ),
#            "dBdx"            => ( d[:∂B∂x].odata,                            ("Nx", "Ny", "Nz_c") ),
#            "dBdy"            => ( d[:∂B∂y].odata,                            ("Nx", "Nyp1", "Nz_c") ),
            #"u_U"             => ( d[:u_U].odata,                             ("Nx", "Ny",   "Nz_f") ),
            #"v_V"             => ( d[:v_V].odata,                             ("Nx", "Nyp1", "Nz_f") ),
#            "w_W"            => ( d[:w_W].odata,                             ("Nx", "Ny", "Nz_fp1") ),
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

