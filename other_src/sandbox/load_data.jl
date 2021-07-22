using NCDatasets
using Formatting
using Statistics

nonan = (x) -> replace(x, missing => NaN)


Dataset("raw_sfc_data.singlepoint.nc") do ds

    global τx_data     = (mean( reshape( ds["TAUX"][1, 1, :], 12, :), dims=(2,))[:, 1] |> nonan) * 0.1   # dyne/cm^2 = 0.1 N/m^2
    global τy_data     = (mean( reshape( ds["TAUY"][1, 1, :], 12, :), dims=(2,))[:, 1] |> nonan) * 0.1
    global swflx_data  = mean( reshape( - ds["SHF_QSW"][1, 1, :], 12, :), dims=(2,))[:, 1] |> nonan
    global nswflx_data = mean( reshape( - (ds["SHF"][1, 1, :] - ds["SHF_QSW"][1, 1, :]), 12, :), dims=(2,))[:, 1] |> nonan


    
end

Dataset("qflx.singlepoint.nc") do ds
    global qflx_data = ds["qdp"][1,1,:]  |> nonan
    global h_ML_data = ds["hblt"][1,1,:] |> nonan
end

