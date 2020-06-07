using NCDatasets
using ArgParse
using JSON
using Statistics
using Formatting

function parse_commandline()

    s = ArgParseSettings()
    @add_arg_table s begin

        "--data-file-prefix"
            help = "Data filename prefix including folder and path until the timestamp. File extension `nc` is assumed."
            arg_type = String
            required = true
 
        "--domain-file"
            help = "Domain file."
            arg_type = String
            default = ""

        "--input-file"
            help = "Input Qflx file name."
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

if doy != 365
    throw(ErrorException("doy is not 365! Now it is " * string(doy)))
end

fday = zeros(Int64, 12)
fday[1] = 1
for m=2:12
    fday[m] = fday[m-1] + dom[m-1]
end

nomissing = (x) -> replace(x, missing=>NaN)

if parsed["domain-file"] == ""
    parsed["domain-file"] = parsed["input-file"]
end

Dataset(parsed["domain-file"], "r") do ds

    global mask, xc, yc, area, ni, nj
    mask = ds["mask"][:] |> nomissing
    xc   = ds["xc"][:]   |> nomissing
    yc   = ds["yc"][:]   |> nomissing
    area = ds["area"][:] |> nomissing
    ni, nj = size(area)



end

idx       = mask .== 1.0
_area     = area[idx]
_area_sum = sum(_area)

function calAreaAvg(field, ifrac=nothing)
    if ifrac == nothing
        return sum(_area .* (field[idx])) / _area_sum
    else
        area_ifrac = _area .* ifrac[idx]
        return sum(area_ifrac .* (field[idx])) / sum(area_ifrac)
    end
end

function calDailyAreaAvg(field, ifrac=nothing)
    if ifrac == nothing
        return convert(Array{Float64}, [calAreaAvg(view(field, :, :, d)) for d=1:doy])
    else
        return convert(Array{Float64}, [calAreaAvg(view(field, :, :, d), ifrac) for d=1:doy])
    end
end

function adjust!(;
    Qflx  :: AbstractArray{Float64, 2}
    ifrac :: AbstractArray{Float64, 2}
)

    aavg_all_Qflx = calAreaAvg(Qflx) 
    
    aavg_openocean_Qflx = calAreaAvg(Qflx, ifrac) 
    @. Qflx -= aavg_openocean_Qflx * ( 1.0 - ifrac )
    
    aavg_lastall_Qflx = calAreaAvg(Qflx) 
    @. Qflx -= aavg_lastall_Qflx

    return aavg_all_Qflx, aavg_openocean_Qflx, aavg_lastall_Qflx
end




Dataset(parsed["input-file"], "r") do ds

    global Qflx_T_old, Qflx_S_old
    
    Qflx_T_old = ds["Qflx_T"][:] |> nomissing
    Qflx_S_old = ds["Qflx_S"][:] |> nomissing

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

    global varmap = Dict(
        "Qflx_T_correction" => "qflx_T_correction",
        "Qflx_S_correction" => "qflx_S_correction",
        "TSAS_clim"         => "TSAS_clim",
        "SSAS_clim"         => "SSAS_clim",
        "ifrac"             => "ifrac",
    )
    global data = Dict()

    for (key, varname) in varmap
        data[key] = zeros(Float64, ni, nj, doy)
    end

    global cnt = zeros(Float64, doy)  # this is necessary

    for m=1:12

        println(format("Doing month : {:02d}", m))

        for y=parsed["beg-year"]:parsed["end-year"]

            tmp_data = Dict()
            Dataset(format("{:s}{:04d}-{:02d}.nc", parsed["data-file-prefix"], y, m)) do ds
                
                for (key, varname) in varmap
                    tmp_data[key] = ds[varname][:] |> nomissing
                end
            end

            _, _, ndays = size(tmp_data["Qflx_T_correction"])

            if y == 1 && m == 1
                println("Skip the first month of the very first year.")
                continue
            end
            
            if ndays != dom[m]
                throw( ErrorException( format("Month {:02d} does not contain expected days. {:d} instead of {:0d} days.", m, ndays, dom[m])))
            end

            beg_day = fday[m]
            end_day = beg_day + dom[m] - 1

            for (key, varname) in varmap
                data[key][:, :,beg_day:end_day] += tmp_data[key]
            end
            
            cnt[beg_day:end_day] .+= 1.0

        end

    end

    if sum(data["Qflx_T_correction"]) == 0
        throw(ErrorException("T!!!"))
    end

    if sum(data["Qflx_S_correction"]) == 0
        throw(ErrorException("S!!!"))
    end

    println("sum of S correction: ", sum(data["Qflx_S_correction"]))

    for m = 1:12
        beg_day = fday[m]
        end_day = beg_day + dom[m] - 1
        println("cnt of month ", m, " = ", cnt[beg_day:end_day])
    end

    if any(cnt .== 0)
        throw( ErrorException("Some of cnt is zero!") )
    end
 
    # Average
    for (key, varname) in varmap
        var = data[key]
        for d=1:doy
            var[:, :, d] ./= cnt[d]

        end
    end    

    # Transform into energy flux and salt flux
    data["TSAS_clim"] .*= 3996.0 * 1026.0

end



# calculate average information
let

    global avg     = Dict()
    
    global aavg_all_Qflx_T = zeros(Float64, doy)
    global aavg_all_Qflx_S = zeros(Float64, doy)
    
    global aavg_openocean_Qflx_T = zeros(Float64, doy)
    global aavg_openocean_Qflx_S = zeros(Float64, doy)
    
    global aavg_lastall_Qflx_T = zeros(Float64, doy)
    global aavg_lastall_Qflx_S = zeros(Float64, doy)

    global Qflx_T_final = Qflx_T_old + data["Qflx_T_correction"]
    global Qflx_S_final = Qflx_S_old + data["Qflx_S_correction"]

    for d=1:doy
         
        aavg_all_Qflx_T[d], aavg_openocean_Qflx_T[d], aavg_lastall_Qflx_T[d] = adjust!(
            Qflx  = view(Qflx_T_final :, :, d),
            ifrac = view(data["ifrac"], :, :, d),
        )

        aavg_all_Qflx_S[d], aavg_openocean_Qflx_S[d], aavg_lastall_Qflx_S[d] = adjust!(
            Qflx  = view(Qflx_S_final :, :, d),
            ifrac = view(data["ifrac"], :, :, d),
        )

    end


    lnd_idx = mask .== 0.0
    for var in [Qflx_T_final, Qflx_S_final]
        for d = 1:doy
            view(var, :, :, d)[lnd_idx] .= NaN
        end
    end

    global aavg_Qflx_T_old = calDailyAreaAvg(Qflx_T_old)
    global aavg_Qflx_S_old = calDailyAreaAvg(Qflx_S_old)

    global aavg_Qflx_T_final = calDailyAreaAvg(Qflx_T_final)
    global aavg_Qflx_S_final = calDailyAreaAvg(Qflx_S_final)

    for d = 1:doy
        println(
            format("[{:03d}] Qflx_T: all_avg={:.5e}, openocean_avg={:.5e}, lastall_avg={:.5e}. Double check final avg: {:.5e}. SAS_clim={.5e} ",
                d,
                aavg_all_Qflx_T[d],
                aavg_openocean_Qflx_T[d],
                aavg_lastall_Qflx_T[d],
                aavg_Qflx_T_final[d],
                avg["TSAS_clim"][d],
            )
        )
    end

    for d = 1:doy
        println(
            format("[{:03d}] Qflx_S: all_avg={:.5e}, openocean_avg={:.5e}, lastall_avg={:.5e}. Double check final avg: {:.5e}. SAS_clim={.5e} ",
                d,
                aavg_all_Qflx_S[d],
                aavg_openocean_Qflx_S[d],
                aavg_lastall_Qflx_S[d],
                aavg_Qflx_S_final[d],
                avg["SSAS_clim"][d],
            )
        )
    end

end

println("========================")
println("mean_Qflx_T_old: ", mean(aavg_Qflx_T_old))
println("mean_Qflx_S_old: ", mean(aavg_Qflx_S_old))
println("========================")
println("mean_Qflx_T_final: ", mean(aavg_Qflx_T_final)
println("mean_Qflx_S_final: ", mean(aavg_Qflx_S_final)
println("========================")


Dataset(parsed["output-file"], "c") do ds

    defDim(ds, "ni", ni)
    defDim(ds, "nj", nj)
    defDim(ds, "time", doy)

    for (varname, vardata, vardim, attrib) in [
        ("TSAS_clim", data["TSAS_clim"], ("ni", "nj", "time"), Dict()),
        ("SSAS_clim", data["SSAS_clim"], ("ni", "nj", "time"), Dict()),
        ("Qflx_T_correction", data["Qflx_T_correction"], ("ni", "nj", "time"), Dict()),
        ("Qflx_S_correction", data["Qflx_S_correction"], ("ni", "nj", "time"), Dict()),
        ("Qflx_T_final", Qflx_T_final, ("ni", "nj", "time"), Dict()),
        ("Qflx_S_final", Qflx_S_final, ("ni", "nj", "time"), Dict()),
        ("area", area, ("ni", "nj"), Dict()),
        ("mask", mask, ("ni", "nj"), Dict()),
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



