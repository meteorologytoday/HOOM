using Statistics
using NCDatasets

using ArgParse
using JSON
using Formatting

include("LinearRegression.jl")
include("nanop.jl")

correlation = (x1, x2) -> x1' * x2 / (sum(x1.^2)*sum(x2.^2)).^0.5

function parse_commandline()

    s = ArgParseSettings()
    @add_arg_table s begin

        "--data-file"
            help = "Ocean data file. New variable will be appended."
            arg_type = String
            required = true
 
        "--domain-file"
            help = "Domain file."
            arg_type = String
            required = true
 
        "--output-file"
            help = "Output file."
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

        "--varname"
            help = "Variable name. If it is 2D then it should be (lon, lat, time), if it is 3D then it should be (lon, lat, lev/depth, time)"
            arg_type = String
            required = true

              
        "--detrend"
            help = "Whether to detrend the data or not."
            action = :store_true 


     
    end

    return parse_args(ARGS, s)
end

parsed = parse_commandline()
print(json(parsed, 4))

output_file = parsed["output-file"]

Dataset(parsed["data-file"], "r") do ds

    global z_exists

    var = ds[parsed["varname"]]
    dims = size(var)

    global beg_t = (parsed["beg-year"] - 1) * 12 + 1
    global end_t = (parsed["end-year"] - 1) * 12 + 12


    if length(dims) == 3 # 2D + time
        z_exists = false
    elseif length(dims) == 4 # 3D + time
        z_exists = true
    else
        ErrorException("Unknown dimension") |> throw
    end




    if z_exists
        global (Nx, Ny, Nz, _ ) = size(var)
    else
        global (Nx, Ny, _ ) = size(var)
        Nz = 1
    end

    global Nt = end_t - beg_t + 1
    
    #println("Selected Range: ", rng)
    #println("Nt: ", Nt)

    
    if mod(Nt, 12) != 0
        ErrorException("Time record is not multiple of 12") |> throw
    end
    
    global nyears = Int64(Nt / 12)
end

Dataset(parsed["domain-file"], "r") do ds
    global mask = replace(ds["mask"], missing=>NaN)
end

# MM   = Monthly Mean
# A    = Anomaly
# YYC  = Year-to-Year Correlation

data_MM    = zeros(Float64, Nx, Ny, Nz, 12)
data_MA    = zeros(Float64, Nx, Ny, Nz, Nt)
data_MAVAR = zeros(Float64, Nx, Ny, Nz, 12)

x = collect(Float64, 1:Nt)
ds = Dataset(parsed["data-file"], "r")
for k=1:Nz
    rng = ( z_exists ) ? (:, :, k, beg_t:end_t) : (:, :, beg_t:end_t)
    global data  = replace(ds[parsed["varname"]][rng...], missing=>NaN)

    for i=1:Nx, j=1:Ny

        d = view(data, i, j, :)
        if parsed["detrend"]
            d = detrend(x, data)
        end

        data_MM[i, j, k, :] = mean( reshape(d, 12, :), dims=2 )[:, 1]
        data_MA[i, j, k, :] = d - repeat( data_MM[i, j, k, :], outer=nyears)
        
        for m = 1:12
            d_yy = view(data_MA, i, j, k, m:12:(m+nyears*12-1))
            data_MAVAR[i, j, k, m] = var(d_yy)
        end
    end     
end

close(ds)


Dataset(output_file, "c") do ds

    defDim(ds, "months", 12)
    defDim(ds, "time", Inf)
    defDim(ds, "Nx", Nx)
    defDim(ds, "Ny", Ny)
    defDim(ds, "Nz", Nz)
    
    for (varname, vardata, vardim, attrib) in [
        (format("{:s}_MM",    parsed["varname"]),       data_MM,    ("Nx", "Ny", "Nz", "months"), Dict()),
    #    (format("{:s}_MA",    parsed["varname"]),      data_MA,    ("Nx", "Ny", "Nz", "time"),   Dict()),
        (format("{:s}_MAVAR", parsed["varname"]),       data_MAVAR, ("Nx", "Ny", "Nz", "months"), Dict()),
        (format("{:s}_ZONAL_MEAN", parsed["varname"]),  nanmean( data_MM,    dims=(1,) )[1, :, :, :], ("Ny", "Nz", "months"), Dict()),
        (format("{:s}_ZONAL_MAVAR", parsed["varname"]), nanmean( data_MAVAR, dims=(1,) )[1, :, :, :], ("Ny", "Nz", "months"), Dict()),
    ]

        if ! haskey(ds, varname)
            var = defVar(ds, varname, Float64, vardim)
            var.attrib["_FillValue"] = 1e20
        end

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
