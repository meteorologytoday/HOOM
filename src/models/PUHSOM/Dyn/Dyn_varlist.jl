function getCoordinateVariable(m::DynModel)
        return Dict(
            "f"               => ( m.env.gi.c_f,                                                ("Nx", "Ny",) ),
        )
end

function getCompleteVariableList(m::DynModel)
        s = m.state
        c = m.core
        return Dict(
            "Phi"               => ( s.Φ,                                                ("Nx", "Ny",) ),
            "U"                 => ( s.U,                                                ("Nx", "Ny",) ),
            "V"                 => ( s.V,                                                ("Nx", "Nyp1",) ),
            "u_c"               => ( s.u_c,                                              ("Nx", "Ny", "Nz_c") ),
            "v_c"               => ( s.v_c,                                              ("Nx", "Nyp1", "Nz_c") ),
            "u"                 => ( s.u,                                                ("Nx", "Ny", "Nz_c") ),
            "v"                 => ( s.v,                                                ("Nx", "Nyp1", "Nz_c") ),
            "b"                 => ( s.u,                                                ("Nx", "Ny", "Nz_c") ),
            "fu"                 => ( c.wksp.cV[2],                                      ("Nx", "Nyp1", "Nz_c") ),
            "fv"                 => ( c.wksp.cU[2],                                      ("Nx", "Ny", "Nz_c") ),




#=
            # RECORD
            "X"                 => ( PermutedDimsArray(s.X, (2, 3, 1, 4)),             ("Nx", "Ny", "Nz_f", "NX") ),
#            "XFLUX_bot"         => ( m.tcr_adv.XFLUX_bot,                                    ("Nx", "Ny", "NX") ),
#            "XFLUX_DEN_x"       => ( PermutedDimsArray(m.tcr_adv.XFLUX_DEN_x, (2, 3, 1, 4)),   ("Nx", "Ny", "Nz_f", "NX") ),
#            "XFLUX_DEN_y"       => ( PermutedDimsArray(m.tcr_adv.XFLUX_DEN_y, (2, 3, 1, 4)),   ("Nx", "Nyp1", "Nz_f", "NX") ),
#            "XFLUX_DEN_z"       => ( PermutedDimsArray(m.tcr_adv.XFLUX_DEN_z, (2, 3, 1, 4)),   ("Nx", "Ny", "Nz_fp1", "NX") ),
#            "div"               => ( PermutedDimsArray(m.tcr_adv.div, (2, 3, 1)),            ("Nx", "Ny", "Nz_f") ),
            "u_f"               => ( PermutedDimsArray(s.u_f, (2, 3, 1)),              ("Nx", "Ny", "Nz_f") ),
            "v_f"               => ( PermutedDimsArray(s.v_f, (2, 3, 1)),              ("Nx", "Nyp1", "Nz_f") ),
            "w_f"               => ( PermutedDimsArray(s.w_f, (2, 3, 1)),              ("Nx", "Ny", "Nz_fp1") ),
            "u_c"               => ( PermutedDimsArray(s.u_c, (2, 3, 1)),              ("Nx", "Ny", "Nz_c") ),
            "v_c"               => ( PermutedDimsArray(s.v_c, (2, 3, 1)),              ("Nx", "Nyp1", "Nz_c") ),
            "Phi"               => ( s.Φ,                                                ("Nx", "Ny",) ),
            "U"               => ( s.U,                                                ("Nx", "Ny",) ),
            "V"               => ( s.V,                                                ("Nx", "Nyp1",) ),
            "u"               => ( PermutedDimsArray(s.u, (2, 3, 1)),              ("Nx", "Ny", "Nz_c") ),
            "v"               => ( PermutedDimsArray(s.v, (2, 3, 1)),              ("Nx", "Nyp1", "Nz_c") ),
            "b"               => ( PermutedDimsArray(s.b, (2, 3, 1)),              ("Nx", "Ny", "Nz_c") ),
            "u_aux"               => ( PermutedDimsArray(c.u_aux, (2, 3, 1)),              ("Nx", "Ny", "Nz_c") ),
            "v_aux"               => ( PermutedDimsArray(c.v_aux, (2, 3, 1)),              ("Nx", "Nyp1", "Nz_c") ),
=#
        )
end

