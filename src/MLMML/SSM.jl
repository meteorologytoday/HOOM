
include("MLMML.jl")


module SSM

using Printf
using Formatting
using ..MLMML

missing_value = 1e20


include("Workspace.jl")
include("OceanColumnCollection.jl")

function stepOceanColumnCollection!(;
    occ   :: OceanColumnCollection,
    Δt    :: Float64,

)
    wksp = occ.wksp
    
    wksp.fric_u .= sqrt.(sqrt.((wksp.taux).^2.0 + (wksp.tauy).^2.0) / MLMML.ρ)

    wksp.hflx   .*= (MLMML.α * MLMML.g / MLMML.ρ / MLMML.c_p)
    wksp.swflx  .*= (MLMML.α * MLMML.g / MLMML.ρ / MLMML.c_p)
    
    for l = 1:occ.N_ocs

        if occ.mask[l] == 0.0
            continue
        end        

        MLMML.stepOceanColumn!(;
            oc     = occ.ocs[l],
            fric_u = wksp.fric_u[l] * (1.0 - wksp.ifrac[l]),
            B0     = wksp.hflx[l],
            J0     = wksp.swflx[l],
            Δt     = Δt,
        )
    end

end

function maskData!(occ::OceanColumnCollection, arr::Array{Float64})
    for i = 1:occ.N_ocs
        if occ.mask[i] == 0.0
            arr[i] = missing_value
        end
    end
end

function getInfo!(;
    occ      :: OceanColumnCollection,
    sst      :: Union{Array{Float64}, Nothing} = nothing,
    mld      :: Union{Array{Float64}, Nothing} = nothing,
    qflx2atm :: Union{Array{Float64}, Nothing} = nothing,
)
    if mld != nothing
        for l = 1:occ.N_ocs
            if occ.mask[l] == 0.0
                continue
            end
            mld[l] = occ.ocs[l].h_ML
        end
    end

    if sst != nothing

        for l = 1:occ.N_ocs

          if occ.mask[l] == 0.0
              continue
          end
          sst[l] = occ.ocs[l].b_ML / (MLMML.α * MLMML.g) + MLMML.T_ref
        end

    end

    if qflx2atm != nothing
        for l = 1:occ.N_ocs
            if occ.mask[l] == 0.0
                continue
            end
            qflx2atm[l] = occ.ocs[l].qflx2atm
        end
    end

end



end
