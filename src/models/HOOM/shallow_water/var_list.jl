function getCompleteVariableList(m::DynModel)
        s = m.state
        return Dict(
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
            "Phi"               => ( s.Φ,                                                ("Nx", "Ny",) ),
        )
end

