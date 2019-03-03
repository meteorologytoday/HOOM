function getInfo!(;
    occ :: OceanColumnCollection,
    sst :: Union{Array{Float64}, Nothing} = nothing,
    mld :: Union{Array{Float64}, Nothing} = nothing,
    idx :: Union{Integer, Nothing} = nothing,
)
    if mld != nothing
        if idx == nothing
            throw(ErrorException("[getInfo!] `mld` is reqested but `idx` is not given."))
        end

        for l = 1:occ.N_ocs
            if occ.mask[l] == 0.0
                continue
            end
            mld[l] = occ.h_ML[l, idx]
        end
    end

    if sst != nothing

        for l = 1:occ.N_ocs

          if occ.mask[l] == 0.0
              continue
          end
          sst[l] = occ.b_ML[l] / (MLM2L.α * MLM2L.g)
        end

    end
end



