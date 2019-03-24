include("tools.jl")
using NCDatasets
using Printf
using Formatting

α = 3e-4
β = 1e-3
g = 9.8
ρ    = 1027.0  # kg / m^3
g    = 9.8     # m / s^2
c_p  = 3985.0  # J / kg / K


dtype = Float32

ifilename_Temp = "RG_ArgoClim_Temp.nc"
ifilename_Psal = "RG_ArgoClim_Psal.nc"
ds = Dataset(ifilename_Temp,"r")
missing_value = convert(dtype, ds["ARGO_TEMPERATURE_ANOMALY"].attrib["_FillValue"])
lon  = ds["LONGITUDE"][:]
lat  = ds["LATITUDE"][:]
pres  = ds["PRESSURE"][:]
time = collect(dtype, 1:length(ds["TIME"]))
close(ds)

println(typeof(missing_value))


function ρ(PSU, T) 
     SA = 35.16504/35.0 * PSU
end
