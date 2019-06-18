using Statistics
using NCDatasets

using ArgParse
using JSON

correlation = (x1, x2) -> x1' * x2 / (sum(x1.^2)*sum(x2.^2)).^0.5

function parse_commandline()

    s = ArgParseSettings()
    @add_arg_table s begin

        "--data-file"
            help = "Atm data file. New variable will be appended."
            arg_type = String
            required = true
 
        "--domain-file"
            help = "Domain file."
            arg_type = String
            required = true
      
    end

    return parse_args(ARGS, s)
end

parsed = parse_commandline()
print(json(parsed, 4))

output_file = joinpath(dirname(parsed["data-file"]), "atm_anomalies.nc")

Dataset(parsed["data-file"], "r") do ds

    dims = size(ds["PSL"])

    if length(dims) == 3 # 2D + time
        rng = (:,) 
    elseif length(dims) == 4 # 2D + time
        rng = (:,:,1,:) 
    else
        ErrorException("Unknown dimension") |> throw
    end

    global PSL  = replace(ds["PSL"][rng...], missing=>NaN)

    global (Nx, Ny, Nt) = size(PSL)
    
    if mod(Nt, 12) != 0

        PSL = PSL[:, :, 1:12*floor(Integer, Nt/12.0)]
        (Nx, Ny, Nt) = size(PSL)
        
        #ErrorException("Time record is not multiple of 12") |> throw
    end
    
    global nyears = Int64(Nt / 12)
end

Dataset(parsed["domain-file"], "r") do ds
    global mask = replace(ds["mask"], missing=>NaN)
end

# MM   = Monthly Mean
# A    = Anomaly
# YYC  = Year-to-Year Correlation

PSLMM   = zeros(Float64, size(PSL)[1:2]..., 12)
PSLA    = zeros(Float64, size(PSL)...)
PSLAYYC = zeros(Float64, size(PSL)[1:2]..., 12)
PSLAVAR = zeros(Float64, size(PSL)[1:2]..., 12)


for i=1:Nx, j=1:Ny

    d = view(PSL, i, j, :)

    PSLMM[i, j, :] = mean( reshape(d, 12, :), dims=2 )[:, 1]
    PSLA[i, j, :]  = d - repeat( PSLMM[i, j, :], outer=nyears)
    
    for m = 1:12
        d_yy = view(PSLA, i, j, m:12:(m+nyears*12-1))
        PSLAYYC[i, j, m] = correlation(d_yy[1:end-1], d_yy[2:end])
        PSLAVAR[i, j, m] = std(d_yy)
    end
 
end



Dataset(output_file, "c") do ds

    defDim(ds, "months", 12)
    defDim(ds, "time", Inf)
    defDim(ds, "Nx", Nx)
    defDim(ds, "Ny", Ny)
    
    for (varname, vardata, vardim, attrib) in [
        ("PSLMM",  PSLMM, ("Nx", "Ny", "months"), Dict()),
        ("PSLA",    PSLA, ("Nx", "Ny", "time"),   Dict()),
        ("PSLAYYC",  PSLAYYC, ("Nx", "Ny", "months"), Dict()),
        ("PSLAVAR",  PSLAVAR, ("Nx", "Ny", "months"), Dict()),
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
