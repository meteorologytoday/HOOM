if MD.configs["enable_short_term_archive"]

    if MD.configs["daily_record"]
        
        if t_flags["new_month"] && substep == 1
                filename = format("{}.ocn.h.daily.{:04d}-{:02d}.nc", MD.casename, t[1], t[2])
                RecordTool.setNewNCFile!(
                    MD.recorders["daily_record"],
                    joinpath(MD.configs["short_term_archive_dir"], filename)
                )
        end

        RecordTool.record!(
            MD.recorders["daily_record"];
            avg_and_output = ( t_flags["new_day"] && substep==1 && t_cnt != 1)
        )

        appendLine(MD.configs["short_term_archive_list"], filename)
    end
    
    if MD.configs["monthly_record"]

        if t_flags["new_year"] && substep == 1
                filename = format("{}.ocn.h.monthly.{:04d}.nc", MD.casename, t[1])
                RecordTool.setNewNCFile!(
                    MD.recorders["monthly_record"],
                    joinpath(MD.configs["short_term_archive_dir"], filename)
                )
        end

        RecordTool.record!(
            MD.recorders["monthly_record"];
            avg_and_output = ( t_flags["new_month"] && substep==1 && t_cnt != 1)
        )


    end

end

