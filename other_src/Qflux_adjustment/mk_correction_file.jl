using NCDatasets
using ArgParse
using JSON
using Statistics

function parse_commandline()

    s = ArgParseSettings()
    @add_arg_table s begin

        "--input-file"
            help = "File "
            arg_type = String
            required = true
        
        "--domain-file"
            help = "Domain file."
            arg_type = String
            required = true

        "--output-ConH-file"
            help = "Output ConH file name."
            arg_type = String
            required = true


        "--output-CycH-file"
            help = "Output CycH file name."
            arg_type = String
            required = true

    end

    return parse_args(ARGS, s)
end

parsed = parse_commandline()
print(json(parsed, 4))


nomissing = (x) -> replace(x, missing=>NaN)

Dataset(parsed["domain-file"], "r") do ds
    global mask, xc, yc, area, ni, nj
    mask = ds["mask"][:] |> nomissing
    xc   = ds["xc"][:]   |> nomissing
    yc   = ds["yc"][:]   |> nomissing
    area = ds["area"][:] |> nomissing
    ni, nj = size(area)
end


