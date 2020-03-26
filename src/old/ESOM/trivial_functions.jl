
function TS2b(T::Float64, S::Float64)
    return g * (α * (T - T_ref) - β * (S - S_ref))
end

function OC_updateB!(
    occ :: OceanColumnCollection,
    i   :: Integer,
    j   :: Integer,
)

    occ.bs[i, j, 1] = TS2b(occ.Ts[i, j, 1], occ.Ss[i, j, 1])
    occ.bs[i, j, 2] = TS2b(occ.Ts[i, j, 2], occ.Ss[i, j, 2])

end

function updateB!(occ::OceanColumnCollection)

    for i=1:occ.Nx, j=1:occ.Ny
        OC_updateB!(occ, i, j)
    end

end
