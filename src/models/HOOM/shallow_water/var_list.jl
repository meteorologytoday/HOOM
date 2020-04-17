function getCompleteVariableList(m::Model)

        return Dict(
            # RECORD
            "X"                 => ( PermutedDimsArray(m.state.X, (2, 3, 1, 4)),             ("Nx", "Ny", "Nz_f", "NX") ),
            "XFLUX_bot"         => ( m.tcr_adv.XFLUX_bot,                                    ("Nx", "Ny", "NX") ),
            "XFLUX_DEN_x"       => ( PermutedDimsArray(m.tcr_adv.XFLUX_DEN_x, (2, 3, 1, 4)),   ("Nx", "Ny", "Nz_f", "NX") ),
            "XFLUX_DEN_y"       => ( PermutedDimsArray(m.tcr_adv.XFLUX_DEN_y, (2, 3, 1, 4)),   ("Nx", "Nyp1", "Nz_f", "NX") ),
            "XFLUX_DEN_z"       => ( PermutedDimsArray(m.tcr_adv.XFLUX_DEN_z, (2, 3, 1, 4)),   ("Nx", "Ny", "Nz_fp1", "NX") ),
            "div"               => ( PermutedDimsArray(m.tcr_adv.div, (2, 3, 1)),            ("Nx", "Ny", "Nz_f") ),
            "u_f"               => ( PermutedDimsArray(m.state.u_f, (2, 3, 1)),              ("Nx", "Ny", "Nz_f") ),
            "v_f"               => ( PermutedDimsArray(m.state.v_f, (2, 3, 1)),              ("Nx", "Nyp1", "Nz_f") ),
            "w_f"               => ( PermutedDimsArray(m.state.w_f, (2, 3, 1)),              ("Nx", "Ny", "Nz_fp1") ),
        )
end

