function getCoordinateVariable(m::DynModel)
        return Dict(
            "f"               => ( m.env.gi.c_f,                         ("Nx", "Ny",) ),
        )
end

function getCompleteVariableList(m::DynModel)
        s = m.state
        c = m.core
        return Dict(
            "Phi"               => ( s.Φ,                                ("Nx", "Ny",) ),
            "U"                 => ( s.U,                                ("Nx", "Ny",) ),
            "V"                 => ( s.V,                                ("Nx", "Nyp1",) ),
            "u_total"           => ( s.u_total,                          ("Nx", "Ny", "Nz_c") ),
            "v_total"           => ( s.v_total,                          ("Nx", "Nyp1", "Nz_c") ),
            "u"                 => ( s.u,                                ("Nx", "Ny", "Nz_c") ),
            "v"                 => ( s.v,                                ("Nx", "Nyp1", "Nz_c") ),
            "B"                 => ( s.B,                                ("Nx", "Ny", "Nz_c") ),
            "fu"                => ( c.wksp.cV[2],                       ("Nx", "Nyp1", "Nz_c") ),
            "fv"                => ( c.wksp.cU[2],                       ("Nx", "Ny", "Nz_c") ),
        )
end

