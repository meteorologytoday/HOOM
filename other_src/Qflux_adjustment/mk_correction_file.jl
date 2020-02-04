using NCDatasets
using ArgParse
using JSON
using Statistics

function parse_commandline()

    s = ArgParseSettings()
    @add_arg_table s begin

        "--data-file-prefix"
            help = "Data filename prefix including folder and path until the timestamp. File extension `nc` is assumed."
            arg_type = String
            required = true
 
        "--input-file"
            help = "File "
            arg_type = String
            required = true
        
        "--domain-file"
            help = "Domain file."
            arg_type = String
            required = true

        "--output-file"
            help = "Output file name."
            arg_type = String
            required = true

        "--beg-year"
            help = "Begin year of the data."
            arg_type = Int64
            required = true

        "--end-year"
            help = "End year of the data."
            arg_type = Int64
            required = true




    end

    return parse_args(ARGS, s)
end

parsed = parse_commandline()
print(json(parsed, 4))

dom = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
doy = sum(dom)

fday = zeros(12)
fday[1] = 1
for m=2:12
    fday[m] = fday[m-1] + dom[m-1]
end

nomissing = (x) -> replace(x, missing=>NaN)

Dataset(parsed["domain-file"], "r") do ds
    global mask, xc, yc, area, ni, nj
    mask = ds["mask"][:] |> nomissing
    xc   = ds["xc"][:]   |> nomissing
    yc   = ds["yc"][:]   |> nomissing
    area = ds["area"][:] |> nomissing
    ni, nj = size(area)
end

let

    #
    # Variable `cnt` is necessary because CESM simulation skips the first two days
    # to initialize the model. A day goes for the startup run for CAM ( details see
    # http://www.cesm.ucar.edu/models/cesm1.1/cesm/doc/usersguide/c1128.html ), 2nd
    # day in DOCN (see code).
    #
    # So days of Jan of the first year of each run would be 29, 30, 30, 30, ... etc.
    #

    global qflx_correction_mean = zeros(Float64, ni, nj, doy)
    global cnt = zeros(Float64, doy)  # this is nece

    for m=1:12
        for y=parsed["beg-year"]:parsed["end-year"]

            Dataset(format("{:s}{:04d}-{:02d}.nc", parsed["data-file-prefix"], y, m)) do ds

                qflx_correction = ds["qflx_correction"][:] |> nomissing
                _, _, ndays = size(qflx_correction)

                if m == 1
                    if ndays != 31 && ( ! ndays in (29, 30) )
                        throw( ErrorException("January does not contains expected days. (Not 29, 30, or 31 days)") )
                    end
                    
                    wedge_size = 31 - ndays
                    beg_day = wedge_size+1
                    end_day = 31
                    qflux_correction_mean[:, :, beg_day:end_day ] += qflx_correction
                    cnt[beg_day:end_day] .+= 1.0
                    
                else
                    if ndays != dom[m]
                        throw( ErrorException( format("Month {:02d} does not contain expected days. {:d} instead of {:0d} days.", m, ndays, dom[m])))
                    end

                    
                    beg_day = fday[m]
                    end_day = beg_day + dom[m] - 1
                    qflux_correction_mean[:, :, beg_day:end_day ] += qflx_correction
                    cnt[beg_day:end_day] .+= 1.0
 
                end

                

            end
        end
    end
 
    for d=1:doy
        if cnt[d] != 0.0
            qflux_correction_mean[:, :, d] ./= cnt[d]
        end
    end

end

Dataset(parsed["output-file"], "c") do ds

    defDim(ds, "ni", ni)
    defDim(ds, "nj", nj)
    defDim(ds, "time", Inf)

    for (varname, vardata, vardim, attrib) in [
        ("qflx_correction_mean", qflx_correction_mean, ("ni", "nj", "time"), Dict()),
        ("cnt", cnt, ("time", ), Dict()),
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



