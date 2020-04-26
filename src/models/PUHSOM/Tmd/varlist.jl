function getCoordinateVariable(m::TmdModel)
        return Dict(
            "f"               => ( m.env.gi.c_f,                                                ("Nx", "Ny",) ),
        )
end

function getCompleteVariableList(m::TmdModel)
    s = m.state
    c = m.core
    return Dict(
        "X"               => ( PermutedDimsArray(s.X, (2, 3, 1, 4)),                ("Nx", "Ny", "Nz", "NX") ),
        "X_ML"            => ( s.X_ML,                                              ("Nx", "Ny",       "NX") ),
        "u_U"             => ( PermutedDimsArray(s.u_U, (2, 3, 1)),                 ("Nx", "Ny", "Nz"      ) ),
    )
end

