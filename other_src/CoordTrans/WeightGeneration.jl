using Distributed

@everywhere module WeightGeneration

    using NCDatasets
    using SharedArrays
    using Formatting
    using Distributed

    mutable struct GridInfo{T <: AbstractFloat}
        
        N      :: Int64

        gc_lon :: AbstractArray{T}
        gc_lat :: AbstractArray{T}

        area   :: AbstractArray{T}
        mask   :: AbstractArray{T}

        unit_of_angle :: Symbol

        dims   :: AbstractArray{Int64} 

        function GridInfo{T}(;
            gc_lon :: AbstractArray{T,1},
            gc_lat :: AbstractArray{T,1},
            area   :: AbstractArray{T,1},
            mask   :: AbstractArray{T,1},
            unit_of_angle :: Symbol,
            dims   = nothing,
        ) where T <: AbstractFloat 

            N = length(gc_lon)

            for var in [gc_lat, area, mask]
                if length(var) != N
                    throw(ErrorException("Not all input has the same length."))
                end
            end

            if unit_of_angle == :deg

                gc_lon .*= π / 180.0
                gc_lat .*= π / 180.0

            elseif unit_of_angle == :rad
                # do nothing

            else
                throw(ErrorException("`unit_of_angle` must be `:deg` or `:rad`."))
            end

            if dims == nothing
                dims = (N,)
            end

            dims = convert(Array{Int64}, dims)

            if reduce(*, dims) != N
                throw(ErrorException("Dims does not match the number of elements."))
            end
            



            return new(N, gc_lon, gc_lat, area, mask, unit_of_angle, dims)
        end


    end

    function genWeight_NearestNeighbors(
        filename :: AbstractString,
        gi_s     :: GridInfo,
        gi_d     :: GridInfo,
        NNN_max  :: Integer;
        missing_value :: T = 1e20,
    ) where T <: AbstractFloat

        
        trans = SharedArray{T}(NNN_max, gi_d.N)

        # s_coord and d_coord are the coordinates of grid points
        # in 3-dimensional cartesian coordinate

        s_coord = SharedArray{T}(3, gi_s.N)
        d_coord = SharedArray{T}(3, gi_d.N)

        s_NaN_idx = (gi_s.mask .== 0)

        @sync @distributed for i = 1:gi_s.N

            s_coord[1, i] = cos(gi_s.gc_lat[i]) * cos(gi_s.gc_lon[i])
            s_coord[2, i] = cos(gi_s.gc_lat[i]) * sin(gi_s.gc_lon[i])
            s_coord[3, i] = sin(gi_s.gc_lat[i])

        end

        @sync @distributed for i = 1:gi_d.N

            d_coord[1, i] = cos(gi_d.gc_lat[i]) * cos(gi_d.gc_lon[i])
            d_coord[2, i] = cos(gi_d.gc_lat[i]) * sin(gi_d.gc_lon[i])
            d_coord[3, i] = sin(gi_d.gc_lat[i])

        end

        #s_NaN_idx = (s_mask .== 0)

        println("Start making transform matrix... ")

        @time @sync @distributed for i = 1:gi_d.N

            # For every point find its nearest-neighbors

            #print("\r", i, "/", d_N)

            if gi_d.mask[i] == 0
                trans[:, i] .= NaN
                continue
            end

            dist2 = (  (s_coord[1, :] .- d_coord[1, i]).^2
                     + (s_coord[2, :] .- d_coord[2, i]).^2
                     + (s_coord[3, :] .- d_coord[3, i]).^2 )


            # Decided not to apply this condition because in 
            # extreme cases there might be a small area of water
            # that is surrounded by lands.

            dist2[s_NaN_idx] .= NaN
         
            idx_arr = collect(1:gi_s.N)
            sort!(idx_arr; by=(k)->dist2[k])
            trans[:, i] = idx_arr[1:NNN_max]

        end

        println(typeof(gi_s.dims))
        Dataset(filename, "c") do ds

            defDim(ds, "s_N", gi_s.N)
            defDim(ds, "d_N", gi_d.N)
            defDim(ds, "NNN_max", NNN_max)
            defDim(ds, "s_dims", length(gi_s.dims))
            defDim(ds, "d_dims", length(gi_d.dims))

            for (varname, vardata, vardims) in (
                ("NN_idx",    trans, ("NNN_max", "d_N")),
                ("s_gc_lat",  gi_s.gc_lat, ("s_N",)),
                ("s_gc_lon",  gi_s.gc_lon, ("s_N",)),
                ("d_gc_lat",  gi_d.gc_lat, ("d_N",)),
                ("d_gc_lon",  gi_d.gc_lon, ("d_N",)),
                ("s_dims",    gi_s.dims, ("s_dims",)),
                ("d_dims",    gi_d.dims, ("d_dims",)),
                ("s_wgt",     gi_s.area, ("s_N",)),
            )

                print(format("Output data: {} ...", varname))

                dtype = eltype(vardata)

                v = defVar(ds, varname, eltype(vardata), vardims)

                if dtype <: AbstractFloat
                    v.attrib["_FillValue"] = missing_value
                end

                v[:] = vardata
                println("done.")
            end

            
            
            
        end

    end    


    function convertData!(
        NN_idx  :: AbstractArray{I, 2},
        s_wgt   :: AbstractArray{T, 2},
        s_data  :: AbstractArray{T, 1},
        d_data  :: AbstractArray{T, 1},
    ) where T <: AbstractFloat where I <: Integer


        for i = 1 : length(d_data)

            d_data[i] = 0
            wgt_sum = 0.0

            for j = 1:NNN

                idx = NN_idx[j, i]

                data = s_data[idx]
        
                if isfinite(data)
                    wgt_sum += s_wgt[idx]
                    d_data[i] += data
                else
                    break
                end
            end

            d_data[i] = (NNN_real == 0) NaN : d_data[i] / wgt_sum

        end

    end


    function readWeightFile(wgt_filename :: AbstractString)

        ds = Dataset(wgt_filename, "r")

        NN_idx = ds["NN_idx"][:]
        replace!(NN_idx, missing=>NaN)

        s_wgt = ds["s_wgt"][:]
        replace!(s_wgt, missing=>NaN)


        return NN_idx, s_wgt

    end



    function convertData2D(
        in_filename  :: AbstractString,
        out_filename :: AbstractString,
        wgt_filename :: AbstractString,
        varnames2D   :: Tuple,
    )
        NN_idx, s_wgt = readWeightFile(wgt_filename)
        
        ds_in  = Dataset(in_filename, "r")
        ds_out = Dataset(out_filename, "o")

        for i = 1:N
        end 

    end

end
