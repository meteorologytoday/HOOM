function getCompleteVariableList(ocn::Ocean)

        return Dict(
            "T"                  => ( toXYZ(ocn.Ts, :zxy),          ("Nx", "Ny", "Nz_bone") ),
            "S"                  => ( toXYZ(ocn.Ss, :zxy),          ("Nx", "Ny", "Nz_bone") ),
            "b"                  => ( toXYZ(ocn.bs, :zxy),          ("Nx", "Ny", "Nz_bone") ),
            "T_ML"               => ( ocn.T_ML,                     ("Nx", "Ny") ),
            "S_ML"               => ( ocn.S_ML,                     ("Nx", "Ny") ),
            "dTdt_ent"           => ( ocn.dTdt_ent,                 ("Nx", "Ny") ),
            "dSdt_ent"           => ( ocn.dSdt_ent,                 ("Nx", "Ny") ),
            "TSAS_clim"          => ( ocn.TSAS_clim,                ("Nx", "Ny") ),
            "SSAS_clim"          => ( ocn.SSAS_clim,                ("Nx", "Ny") ),
            "TFLUX_bot"          => ( ocn.TFLUX_bot,                ("Nx", "Ny") ),
            "SFLUX_bot"          => ( ocn.SFLUX_bot,                ("Nx", "Ny") ),
            "SFLUX_top"          => ( ocn.SFLUX_top,                ("Nx", "Ny") ),
            "TFLUX_DIV_implied"  => ( ocn.TFLUX_DIV_implied,        ("Nx", "Ny") ),
            "SFLUX_DIV_implied"  => ( ocn.SFLUX_DIV_implied,        ("Nx", "Ny") ),
            "qflx2atm"           => ( ocn.qflx2atm,                 ("Nx", "Ny") ),
            "qflx2atm_pos"       => ( ocn.qflx2atm_pos,             ("Nx", "Ny") ),
            "qflx2atm_neg"       => ( ocn.qflx2atm_neg,             ("Nx", "Ny") ),
            "h_ML"               => ( ocn.h_ML,                     ("Nx", "Ny") ),
            "h_MO"               => ( ocn.h_MO,                     ("Nx", "Ny") ),
            "nswflx"             => ( ocn.in_flds.nswflx,           ("Nx", "Ny") ),
            "swflx"              => ( ocn.in_flds.swflx,            ("Nx", "Ny") ),
            "frwflx"             => ( ocn.in_flds.frwflx,           ("Nx", "Ny") ),
            "vsflx"              => ( ocn.in_flds.vsflx,            ("Nx", "Ny") ),
            "qflx_T"             => ( ocn.in_flds.qflx_T,           ("Nx", "Ny") ),
            "qflx_S"             => ( ocn.in_flds.qflx_S,           ("Nx", "Ny") ),
            "qflx_T_correction"  => ( ocn.qflx_T_correction,        ("Nx", "Ny") ),
            "qflx_S_correction"  => ( ocn.qflx_S_correction,        ("Nx", "Ny") ),
            "Tclim_feeded"       => ( ocn.in_flds.Tclim,            ("Nx", "Ny") ),
            "Sclim_feeded"       => ( ocn.in_flds.Sclim,            ("Nx", "Ny") ),
            "TEMP"               => ( ocn.TEMP,                     ("Nx", "Ny") ),
            "dTEMPdt"            => ( ocn.dTEMPdt,                  ("Nx", "Ny") ),
            "SALT"               => ( ocn.SALT,                     ("Nx", "Ny") ),
            "dSALTdt"            => ( ocn.dSALTdt,                  ("Nx", "Ny") ),
            "fric_u"             => ( ocn.fric_u,                   ("Nx", "Ny") ),
            "taux"               => ( ocn.τx,                       ("Nx", "Ny") ),
            "tauy"               => ( ocn.τy,                       ("Nx", "Ny") ),
            "TFLUX_DEN_z"        => ( toXYZ(ocn.TFLUX_DEN_z, :zxy), ("Nx", "Ny", "zs_bone") ),
            "SFLUX_DEN_z"        => ( toXYZ(ocn.SFLUX_DEN_z, :zxy), ("Nx", "Ny", "zs_bone") ),
            "div"                => ( toXYZ(ocn.div, :zxy),         ("Nx", "Ny", "Nz_bone") ),
            "w_bnd"              => ( toXYZ(ocn.w_bnd, :zxy),       ("Nx", "Ny", "zs_bone") ),
            "u"                  => ( toXYZ(ocn.u, :zxy),           ("Nx", "Ny", "Nz_bone") ),
            "v"                  => ( toXYZ(ocn.v, :zxy),           ("Nx", "Ny", "Nz_bone") ),
            "TFLUX_CONV"         => ( toXYZ(ocn.TFLUX_CONV, :zxy),  ("Nx", "Ny", "Nz_bone") ),
            "SFLUX_CONV"         => ( toXYZ(ocn.SFLUX_CONV, :zxy),  ("Nx", "Ny", "Nz_bone") ),

            #################################################################################
            "area"               => ( ocn.mi.area,                  ("Nx", "Ny") ),
            "mask"               => ( ocn.mi.mask,                  ("Nx", "Ny") ),
            "frac"               => ( ocn.mi.frac,                  ("Nx", "Ny") ),
            "c_lon"              => ( ocn.mi.xc,                    ("Nx", "Ny") ),
            "c_lat"              => ( ocn.mi.yc,                    ("Nx", "Ny") ),
            "zs_bone"            => ( ocn.zs_bone,                  ("zs_bone",) ),
            "topo"               => ( ocn.topo,                     ("Nx", "Ny") ),
            "fs"                 => ( ocn.fs,                       ("Nx", "Ny") ),
            "epsilons"           => ( ocn.epsilons,                 ("Nx", "Ny") ),

        )
end

function getVarDesc(varname)
    return haskey(HOOM.var_desc, varname) ? HOOM.var_desc[varname] : Dict()
end

function getVariableList(ocn::Ocean, keyword::Symbol)

        all_varlist = getCompleteVariableList(ocn)

        output_varnames = []


        if keyword == :ALL

            append!(output_varnames, keys(all_varlist))

        elseif keyword == :ESSENTIAL
            
            append!(output_varnames, [
                "T", "S", "T_ML", "S_ML",
                "TSAS_clim", "SSAS_clim",
                "TFLUX_bot", "SFLUX_bot",
                "SFLUX_top",
                "TFLUX_DIV_implied", "SFLUX_DIV_implied",
                "qflx2atm_pos", "qflx2atm_neg",
                "h_ML",
                "nswflx", "swflx", "frwflx", "vsflx",
                "qflx_T", "qflx_S",
                "qflx_T_correction", "qflx_S_correction",
                "Tclim_feeded", "Sclim_feeded",
                "TEMP", "dTEMPdt", "SALT", "dSALTdt",
                "fric_u", "taux", "tauy",
                "TFLUX_DEN_z", "SFLUX_DEN_z",
                "div",
                "w_bnd", "u", "v", "TFLUX_CONV", "SFLUX_CONV",
            ])

        elseif keyword == :COORDINATE

            append!(output_varnames, [
                "area",
                "mask",
                "frac"
                "c_lon",
                "c_lat",
                "zs_bone", 
            ])

        elseif keyword == :SNAPSHOT_COORDINATE

            append!(output_varnames, [
                "area",
                "mask",
                "frac"
                "c_lon",
                "c_lat",
                "zs_bone", 
                "epsilons",
                "fs",
                "topo",
            ])

        elseif keyword == :SNAPSHOT_RECORD
            
            # These are variables used in snapshot in order
            # to be restored for restart run Record.
            
            append!(output_varnames, [
                "T",
                "S",
                "b",
                "T_ML",
                "S_ML",
                "dTdt_ent",
                "dSdt_ent",
                "TSAS_clim",
                "SSAS_clim",
                "TFLUX_bot",
                "SFLUX_bot",
                "SFLUX_top",
                "TFLUX_DIV_implied",
                "SFLUX_DIV_implied",
                "qflx2atm",
                "qflx2atm_pos",
                "qflx2atm_neg",
                "h_ML",
                "h_MO",
                "nswflx",
                "swflx",
                "frwflx",
                "vsflx",
                "qflx_T",
                "qflx_S",
                "qflx_T_correction",
                "qflx_S_correction",
                "Tclim_feeded",
                "Sclim_feeded",
                "TEMP",
                "dTEMPdt",
                "SALT",
                "dSALTdt",
                "fric_u",
                "taux",
                "tauy",
                "TFLUX_DEN_z",
                "SFLUX_DEN_z",
                "div",
                "w_bnd",
                "u",
                "v",
                "TFLUX_CONV",
                "SFLUX_CONV",
            ])
        end

        output_varlist = Dict()
        for varname in output_varnames
            output_varlist[varname] = all_varlist[varname]
        end

        return output_varlist
end



