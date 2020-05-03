
"""

  This function only sets the value in mixed layer without
  conserving quantity 

"""
 function setMixedLayer!(;
    Ts   :: AbstractArray{Float64, 1},
    Ss   :: AbstractArray{Float64, 1},
    zs   :: AbstractArray{Float64, 1},
    Nz   :: Integer,
    T_ML :: Float64,
    S_ML :: Float64,
    h_ML :: Float64,
)
    FLDO = getFLDO(zs=zs, h_ML=h_ML, Nz=Nz)

    if FLDO > 1
        Ts[1:FLDO-1] .= T_ML
        Ss[1:FLDO-1] .= S_ML
    elseif FLDO == -1
        Ts[1:Nz] .= T_ML
        Ss[1:Nz] .= S_ML
    end
   
    return FLDO 
end


function OC_setMixedLayer!(
    m    :: TmdModel,
    i    :: Integer,
    j    :: Integer;
    T_ML :: Float64,
    S_ML :: Float64,
    h_ML :: Float64,
)

    m.state.h_ML[i, j] = h_ML
    m.state.T_ML[i, j] = T_ML
    m.state.S_ML[i, j] = S_ML
    m.state.FLDO[i, j] = setMixedLayer!(
        Ts   = m.core.cols.T[i, j],
        Ss   = m.core.cols.S[i, j],
        zs   = m.core.cols.z_bnd_av[i, j],
        Nz   = m.env.Nz_av[i, j],
        T_ML = T_ML,
        S_ML = S_ML,
        h_ML = h_ML,
    )

end


function OC_setBuoyancy!(
    m    :: TmdModel,
    i    :: Integer,
    j    :: Integer;
    b    :: AbstractArray{Float64,1},
    b_ML :: Float64,
    h_ML :: Float64,
)

    m.state.b[i, j, :] = b
    OC_setMixedLayer!(m, i, j; b_ML=b_ML, h_ML=h_ML)

end


