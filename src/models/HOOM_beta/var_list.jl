function getCompleteVariableList(mb::ModelBlock)

        return Dict(

            # RECORD
            "TEMP"               => ( mb.fi.sv[:TEMP], :T ),
            "SALT"               => ( mb.fi.sv[:SALT], :T ),
            "UVEL"               => ( mb.fi.sv[:UVEL], :U ),
            "VVEL"               => ( mb.fi.sv[:VVEL], :V ),
            "WVEL"               => ( mb.fi.sv[:WVEL], :W ),
            "TAUX"               => ( mb.fi.TAUX,        :sT ),
            "TAUY"               => ( mb.fi.TAUY,        :sT ),
            "TAUX_east"          => ( mb.fi.TAUX_east,   :sT ),
            "TAUY_north"         => ( mb.fi.TAUY_north,  :sT ),
            "SWFLX"              => ( mb.fi.SWFLX,     :sT ),
            "NSWFLX"             => ( mb.fi.NSWFLX,    :sT ),
            "CHKTEMP"           => ( mb.fi.sv[:CHKTEMP],  :sT ),
            "CHKSALT"           => ( mb.fi.sv[:CHKSALT],  :sT ),
            "ADVT"               => ( mb.fi.sv[:ADVT],   :T ),

            # COORDINATEi
#=
            "area"               => ( ocn.mi.area,                  ("Nx", "Ny") ),
            "mask"               => ( ocn.mi.mask,                  ("Nx", "Ny") ),
            "frac"               => ( ocn.mi.frac,                  ("Nx", "Ny") ),
            "c_lon"              => ( ocn.mi.xc,                    ("Nx", "Ny") ),
            "c_lat"              => ( ocn.mi.yc,                    ("Nx", "Ny") ),
            "zs_bone"            => ( ocn.zs_bone,                  ("NP_zs_bone",) ),
=#
        )
end

function getVarDesc(varname)
    return haskey(HOOM.var_desc, varname) ? HOOM.var_desc[varname] : Dict()
end

function getVariableList(mb::ModelBlock, keywords...)

        all_varlist = getCompleteVariableList(mb)

        output_varnames = Dict()

        function makeSubset(dict, keys)
            subset_dict = Dict()
            for k in keys
                subset_dict[k] = dict[k]
            end
            return subset_dict
        end

        for keyword in keywords
            if keyword == :ALL

                return makeSubset(all_varlist, keys(all_varlist))

            elseif keyword == :QFLX_FINDING

                append!(output_varnames, [
                    "qflx_T", "qflx_S",
                    "qflx_T_correction", "qflx_S_correction",
                    "Tclim", "Sclim", "IFRACclim",
                    "TSAS_clim", "SSAS_clim",
                    "ifrac",
                ])

            elseif keyword == :ESSENTIAL
                
                append!(output_varnames, [
                    "T_ML", "S_ML",
                    "Ts_mixed", "Ss_mixed",
                    "TSAS_clim", "SSAS_clim",
                    "TFLUX_bot", "SFLUX_bot",
                    "SFLUX_top",
                    "TFLUX_DIV_implied", "SFLUX_DIV_implied",
                    "qflx2atm_pos", "qflx2atm_neg",
                    "h_ML", "ifrac",
                    "nswflx", "swflx", "frwflx", "vsflx",
                    "qflx_T", "qflx_S",
                    "seaice_nudge_energy",
                    "TEMP", "dTEMPdt", "SALT", "dSALTdt",
                    "fric_u", "taux", "tauy",
                    "TFLUX_DEN_z", "SFLUX_DEN_z",
                    "div",
                    "w_bnd", "u", "v", "TFLUX_CONV", "SFLUX_CONV",
                    "total_heat", "total_heat_budget_residue",
                    "total_salt", "total_salt_budget_residue",
                ])

            elseif keyword == :DEBUG
                
                append!(output_varnames, [
                    "Ts", "Ss", "T_ML", "S_ML", "b_ML", "bs", "FLDO", "h_MO",
                    "Ts_mixed", "Ss_mixed",
                    "TSAS_clim", "SSAS_clim",
                    "TFLUX_bot", "SFLUX_bot",
                    "SFLUX_top",
                    "TFLUX_DIV_implied", "SFLUX_DIV_implied",
                    "qflx2atm_pos", "qflx2atm_neg",
                    "h_ML", "ifrac",
                    "nswflx", "swflx", "frwflx", "vsflx",
                    "qflx_T", "qflx_S",
                    "seaice_nudge_energy",
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
                    "frac",
                    "c_lon",
                    "c_lat",
                    "zs_bone", 
                ])

            elseif keyword == :BACKGROUND

                append!(output_varnames, [
                    "Ts_clim",
                    "Ss_clim",
                    "h_ML_min",
                    "h_ML_max",
                    "topo",
                    "fs",
                    "epsilons",
                    "K_v",
                    "Dh_T",
                    "Dv_T",
                    "Dh_S",
                    "Dv_S",
                    "we_max",
                    "R",
                    "zeta",
                    "Ts_clim_relax_time",
                    "Ss_clim_relax_time",
                ])

            elseif keyword == :RECORD
                
                # These are variables used in snapshot in order
                # to be restored for restart run Record.
                
                append!(output_varnames, [
                    "Ts",
                    "Ss",
                    "bs",
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
                    "Tclim",
                    "Sclim",
                    "IFRACclim",
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
            else

                throw(ErrorException("Unknown keyword: " * string(keyword)))

            end
        end
        output_varlist = Dict()
        for varname in output_varnames
            output_varlist[varname] = all_varlist[varname]
        end

        return output_varlist
end



