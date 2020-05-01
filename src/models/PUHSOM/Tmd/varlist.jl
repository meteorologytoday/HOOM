function getCoordinateVariable(m::TmdModel)
        return Dict(
            "f"               => ( m.env.gi.c_f,                                                ("Nx", "Ny",) ),
        )
end

function getCompleteVariableList(m::TmdModel)

    @fast_extract m
    
    return Dict(
        "X"               => ( PermutedDimsArray(st.X, (2, 3, 1, 4)),                ("Nx", "Ny", "Nz", "NX") ),
        "X_ML"            => ( st.X_ML,                                              ("Nx", "Ny",       "NX") ),
        "swflx"           => ( fr.swflx,                                             ("Nx", "Ny") ),
        "u_U"             => ( PermutedDimsArray(fr.u_U, (2, 3, 1)),                 ("Nx", "Ny", "Nz"      ) ),
    )
end

