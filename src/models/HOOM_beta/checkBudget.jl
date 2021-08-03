function checkBudget!(
    mb :: ModelBlock,
    Δt :: Float64;
    stage :: Symbol
)


    fi = mb.fi
    co = mb.co

    if stage == :BEFORE_STEPPING

        for i = 1:2
            fi._TMP_CHKX_[:, i] = sum(reshape(co.amo.T_Δz_T * view(fi._X_, :, i), co.gd.Nz, :), dims=1)[1, :]
        end

    elseif stage == :AFTER_STEPPING

        for i = 1:2
            fi._CHKX_[:, i] = sum(reshape(co.amo.T_Δz_T * view(fi._X_, :, i), co.gd.Nz, :), dims=1)[1, :]
        end

        @. fi._CHKX_ = ( fi._CHKX_ - fi._TMP_CHKX_ ) / Δt
     
        fi._CHKX_[:, 1] .-= - ( view(fi.SWFLX, :) + view(fi.NSWFLX, :)) ./ ρcp_sw

        for i = 1:2
            fi._CHKX_[:, i] .-= sum(reshape(co.amo.T_Δz_T * view(fi._ADVX_, :, i), co.gd.Nz, :), dims=1)[1, :]
        end
        
    end      

end
