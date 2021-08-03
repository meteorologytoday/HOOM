function setupForcing!(
    mb :: ModelBlock,
)
    fi = mb.fi
    co = mb.co

    PolelikeCoordinate.project!(
        co.gd, 
        fi.TAUX_east,
        fi.TAUY_north,
        fi.TAUX,
        fi.TAUY;
        direction=:Forward,
        grid=:T,
    )
    
end
