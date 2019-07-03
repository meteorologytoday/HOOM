include("CoordTrans.jl")

using .CoordTrans

using NCDatasets
using Distributed
using SharedArrays
using Formatting
using ArgParse
using JSON

function parse_commandline()

    s = ArgParseSettings()
    @add_arg_table s begin

        "--w-file"
            help = "Generated grid file."
            arg_type = String
            required = true

        "--s-file"
            help = "Source grid file."
            arg_type = String
            required = true
       
        "--d-file"
            help = "Destination grid file."
            arg_type = String
            required = true

        "--s-lat-varname"
            help = "Variable name of latitude in source file."
            arg_type = String
            default = "yc"

        "--s-lon-varname"
            help = "Variable name of longitude in source file."
            arg_type = String
            default = "xc"

        "--s-mask-varname"
            help = "Variable name of mask in source file."
            arg_type = String
            default = "mask"

        "--s-mask-value"
            help = "Mask value representing wanted grids in source file."
            arg_type = String
            required = true


        "--d-lat-varname"
            help = "Variable name of latitude in destination file."
            arg_type = String
            default = "yc"

        "--d-lon-varname"
            help = "Variable name of longitude in destination file."
            arg_type = String
            default = "xc"

        "--d-mask-varname"
            help = "Variable name of mask in destination file."
            arg_type = String
            default = "mask"

        "--d-mask-value"
            help = "Mask value representing wanted grids in destination file."
            arg_type = String
            required = true



    end

    return parse_args(ARGS, s)
end

parsed = parse_commandline()
print(json(parsed, 4))


const NNN_max = 9
const missing_value = 1e20

Dataset(parsed["s-file"], "r") do ds

        global s_mask, s_lat, s_lon

        dims = size(ds[parsed["s-lon-varname"]]) |> collect

        s_mask = replace(reshape(ds[parsed["s-mask-varname"]][:], :), missing=>NaN)

        wanted   = ( s_mask .== parse(Float64, parsed["s-mask-value"]) )
        unwanted = ( s_mask .!= parse(Float64, parsed["s-mask-value"]) )

        s_mask[wanted] .= 1.0
        s_mask[unwanted] .= 0.0

        s_lon  = replace(reshape(ds[parsed["s-lon-varname"]][:], :), missing=>NaN)
        s_lat  = replace(reshape(ds[parsed["s-lat-varname"]][:], :), missing=>NaN)

        global gi_s = CoordTrans.GridInfo(
            gc_lon = s_lon,
            gc_lat = s_lat,
            area   = copy(s_lon),
            mask   = s_mask,
            unit_of_angle = :deg,
            dims = dims,
        )

end

Dataset(parsed["d-file"], "r") do ds
        global d_mask, d_lat, d_lon

        dims = size(ds[parsed["d-lon-varname"]]) |> collect


        d_mask = replace(reshape(ds[parsed["d-mask-varname"]][:], :), missing=>NaN)

        wanted   = ( d_mask .== parse(Float64, parsed["d-mask-value"]) )
        unwanted = ( d_mask .!= parse(Float64, parsed["d-mask-value"]) )

        d_mask[wanted] .= 1.0
        d_mask[unwanted] .= 0.0

        d_lon  = replace(reshape(ds[parsed["d-lon-varname"]][:], :), missing=>NaN)
        d_lat  = replace(reshape(ds[parsed["d-lat-varname"]][:], :), missing=>NaN)

        global gi_d = CoordTrans.GridInfo(
            gc_lon = d_lon,
            gc_lat = d_lat,
            area   = copy(d_lon),  # useless
            mask   = d_mask,
            unit_of_angle = :deg,
            dims = dims,
        )

        d_mask = 1 .- replace(reshape(ds["mask"][:], :), missing=>NaN)
        d_Nx = ds.dim["ni"]
        d_Ny = ds.dim["nj"]
        d_N = d_Nx * d_Ny
        d_lon = replace(reshape(ds["xc"][:], :), missing=>NaN)
        d_lat = replace(reshape(ds["yc"][:], :), missing=>NaN)

end

CoordTrans.genWeight_NearestNeighbors(parsed["w-file"], gi_s, gi_d, NNN_max)

