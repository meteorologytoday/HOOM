
function getConfigDescriptor()

    return Dict(


        :MODEL_MISC => [

            ConfigEntry(
                :init_file,
                :optional,
                [nothing, String],
                nothing
            ),

            ConfigEntry(
                :rpointer_file,
                :optional,
                [String,],
                "";
                desc="If read_restart of the init function is set true, then this entry has to contain a valid rpointer filename.",
            ),

            ConfigEntry(
                :enable_archive,
                :required,
                [Bool,],
            ),


            ConfigEntry(
                :daily_record,
                :optional,
                [AbstractArray, Symbol],
                [],
            ),

            ConfigEntry(
                :monthly_record,
                :optional,
                [AbstractArray, Symbol],
                [],
            ),

        ],
        
        :MODEL_CORE => [

            ConfigEntry(
                :domain_file,
                :required,
                [String,],
            ),

            ConfigEntry(
                :cdata_file,
                :optional,
                [String,],
                nothing,
            ),

            ConfigEntry(
                :cdata_timetype,
                :optional,
                [Any,],
            ),

            ConfigEntry(
                :cdata_beg_time,
                :optional,
                [Any,],
            ),

            ConfigEntry(
                :cdata_end_time,
                :optional,
                [Any,],
            ),

            ConfigEntry(
                :cdata_align_time,
                :optional,
                [Any,],
            ),

            ConfigEntry(
                :z_w,
                :optional,
                [AbstractArray{Float64, 1}],
                nothing;
                desc = "Will be overwritten if :init_file is used.",
            ),

            ConfigEntry(
                :substeps,
                :optional,
                [Integer,],
                8;
            ),

            ConfigEntry(
                :advection_scheme,
                :required,
                [:static, :ekman_codron2012_partition],
            ),

            ConfigEntry(
                :MLD_scheme,
                :required,
                [:prognostic, :datastream, :static],
            ),

            ConfigEntry(
                :Qflx,
                :optional,
                [:on, :off],
                :off,
            ),
            
            ConfigEntry(
                :Qflx_finding,
                :optional,
                [:on, :off],
                :off,
            ),
            
            ConfigEntry(
                :weak_restoring,
                :optional,
                [:on, :off],
                :off,
            ),
            
            ConfigEntry(
                :convective_adjustment,
                :optional,
                [:on, :off],
                :off,
            ),


            ConfigEntry(
                :Ks_H,
                :optional,
                [Float64,],
                1e3;
                desc = "Horizontal tracer diffusivity. Will be overwritten if :init_file is used.",
            ),

            ConfigEntry(
                :Ks_V,
                :optional,
                [Float64,],
                1e-4;
                desc = "Vertical tracer diffusivity. Will be overwritten if :init_file is used.",
            ),

            ConfigEntry(
                :Ks_V_cva,
                :optional,
                [Float64,],
                1.0;
                desc = "Convective adjustment diffusivity. Will be overwritten if :init_file is used.",
            ),

            ConfigEntry(
                :τwk_TEMP,
                :optional,
                [Float64,],
                Inf;
                desc = "Timescale of weak-restoring of temperature. Will be overwritten if :init_file is used",
            ),

            ConfigEntry(
                :τwk_SALT,
                :optional,
                [Float64,],
                Inf;
                desc = "Timescale of weak-restoring of salinity. Will be overwritten if :init_file is used",
            ),

            ConfigEntry(
                :rad_R,
                :optional,
                [Float64,],
                0.58;
                desc = "Fast absorption portion of sunlight as described in Paulson & Simpson (1977). Will be overwritten if :init_file is used",
            ),

            ConfigEntry(
                :rad_ζ1,
                :optional,
                [Float64,],
                0.15;
                desc = "Light penetration length scale as described in Paulson & Simpson (1977). Will be overwritten if :init_file is used",
            ),

            ConfigEntry(
                :rad_ζ2,
                :optional,
                [Float64,],
                23.0;
                desc = "Light penetration length scale as described in Paulson & Simpson (1977). Will be overwritten if :init_file is used",
            ),

            ConfigEntry(
                :ϵ,
                :optional,
                [Float64,],
                1.0 / 86400.0;
            ),

            ConfigEntry(
                :Ekman_layers,
                :optional,
                [Integer,],
                5;
            ),

            ConfigEntry(
                :Returnflow_layers,
                :optional,
                [Integer,],
                5;
            ),
        ],
    )
end