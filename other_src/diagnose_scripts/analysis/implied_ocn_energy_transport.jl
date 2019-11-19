using NCDatasets
using Formatting
using ArgParse
using Statistics
using JSON

include("constants.jl")

function mreplace(x)
    return replace(x, missing=>NaN)
end

function integrate(x, dydx)
    y = copy(dydx) * 0.0

    for i = 2:length(x)
        y[i] = y[i-1] + (dydx[i-1] + dydx[i]) * (x[i] - x[i-1]) / 2.0
    end

    return y
end

function parse_commandline()

    s = ArgParseSettings()
    @add_arg_table s begin

        "--data-file"
            help = "Ocean data file."
            arg_type = String
            required = true
 
        "--output-file"
            help = "Output file."
            arg_type = String
            required = true

        "--domain-file"
            help = "Domain file."
            arg_type = String
            required = true
 
        "--beg-year"
            help = "Year of begin."
            arg_type = Int64
            required = true

        "--end-year"
            help = "Year of end."
            arg_type = Int64
            required = true

    end

    return parse_args(ARGS, s)
end

parsed = parse_commandline()
print(json(parsed, 4))

Dataset(parsed["domain-file"], "r") do ds

    global mask = ds["mask"][:] |> mreplace
    global lat = ds["yc"][1, :] |> mreplace

    global lat_weight = 2π * Re^2.0 * cos.(deg2rad.(lat))
    global ϕ          = deg2rad.(lat)

    global O_idx = (mask .== 0)
    global L_idx = (mask .== 1)

end

Dataset(parsed["data-file"], "r") do ds
 
    beg_t = (parsed["beg-year"] - 1) * 12 + 1
    end_t = (parsed["end-year"] - 1) * 12 + 12
    rng = (:,:,beg_t:end_t) 

    getData = (varname) -> mean( replace(ds[varname][rng...], missing=>0.0, NaN=>0.0), dims=(1, ))[1, :, :]

    global ECONV = getData("TFLUX_DIV_implied") * ρc

    global (Ny, Nt) = size(ECONV)

end

IET_OCN       = zeros(Float64, Ny, Nt)

for t = 1:Nt
    IET_OCN[:, t] = integrate(ϕ, ECONV[:, t] .* lat_weight)
end

Dataset(parsed["output-file"], "c") do ds

    defDim(ds, "Nx",  1)
    defDim(ds, "Ny", Ny)
    defDim(ds, "Nz",  1)
    defDim(ds, "time", Inf)

    for (varname, vardata, vardim, attrib) in [
        ("IET_OCN",     reshape(IET_OCN, 1, Ny, 1, Nt), ("Nx", "Ny", "Nz", "time"), Dict()),
        ("ECONV",       reshape(ECONV, 1, Ny, 1, Nt), ("Nx", "Ny", "Nz", "time"), Dict()),
    ]

        println("Doing var: ", varname)

        var = defVar(ds, varname, Float64, vardim)
        var.attrib["_FillValue"] = 1e20

        var = ds[varname]
        
        for (k, v) in attrib
            var.attrib[k] = v
        end

        rng = []
        for i in 1:length(vardim)-1
            push!(rng, Colon())
        end
        push!(rng, 1:size(vardata)[end])
        var[rng...] = vardata

    end
   

end

