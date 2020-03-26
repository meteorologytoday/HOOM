function makeBlankOceanColumn(;zs::Array{Float64, 1})
    N     = length(zs) - 1
    bs    = zeros(Float64, N)
    K     = 0.0
    b_ML  = 0.0
    h_ML  = h_ML_min
    FLDO  = 1

    oc = OceanColumn(N=N, zs=zs, bs=bs, K=K, b_ML=b_ML, h_ML=h_ML, FLDO=FLDO)
    OC_updateFLDO!(oc)

    return oc
end

function makeSimpleOceanColumn(;
    zs      :: Array{Float64, 1},
    b_slope :: Float64 = 30.0 / 5000.0 * g * α,
    b_ML    :: Float64 = 1.0,
    h_ML    :: Float64 = h_ML_min,
    Δb      :: Float64 = 0.0,
    K       :: Float64 = 1e-5
)

oc = makeBlankOceanColumn(zs=zs)

bs = zeros(Float64, length(zs)-1)
for i = 1:length(bs)
    z = (zs[i] + zs[i+1]) / 2.0
    if z > -h_ML
        bs[i] = b_ML
    else
        bs[i] = b_ML - Δb - b_slope * (-z - h_ML)
    end
end

OC_setBuoyancy!(oc, bs=bs, b_ML=b_ML, h_ML=h_ML)
oc.K = K
OC_updateFLDO!(oc)

return oc
end
