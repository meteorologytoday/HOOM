#!/bin/bash

mkdir -p tmp

set -x

cp NKMLM-gamma_LENS.g37_c2_s500_w100.nc   ./tmp/1.nc
cp NKMLM_LENS.g37_c2_s500_w100.nc         ./tmp/2.nc
cp HMXL_clim_meters.nc                    ./tmp/fom.nc
cp pop_frc.gx3v7.110128.nc                ./tmp/frc.nc
cp docn_forcing.SOM_fixed_MLD_LENS.g37.nc ./tmp/frc.nc

cd tmp

ncks -3 fom.nc fom_v3.nc
ncrename -d nlat,lat -d nlon,lon fom_v3.nc
ncks -O -4 fom_v3.nc fom.nc

ncks -3 frc.nc frc_v3.nc
ncrename -d nj,lat -d ni,lon frc_v3.nc
ncks -O -4 frc_v3.nc frc.nc

ncrename -v Q_mean,Qflux -v h_mean,MLD    1.nc
ncrename -v Q_mean,Qflux -v h_mean,MLD    2.nc
ncrename                 -v HMXL,MLD    fom.nc
ncrename -v qdp,Qflux    -v hblt,MLD    frc.nc

ncap2 -O -s "MLD=MLD;Qflux=-Qflux" frc.nc frc.nc            # qdp is negative of Qflux we calculated

ncbo -O --op_typ=- 1.nc 2.nc   -v Qflux,MLD result_1m2.nc

ncbo -O --op_typ=- 1.nc frc.nc -v Qflux,MLD result_1mfrc.nc
ncbo -O --op_typ=- 2.nc frc.nc -v Qflux,MLD result_2mfrc.nc

ncbo -O --op_typ=- 1.nc fom.nc -v MLD result_1mfom.nc
ncbo -O --op_typ=- 2.nc fom.nc -v MLD result_2mfom.nc

ncrename -v Qflux,Qflux_1m2 -v MLD,MLD_1m2 result_1m2.nc

ncrename -v Qflux,Qflux_1mfrc -v MLD,MLD_1mfrc result_1mfrc.nc
ncrename -v Qflux,Qflux_2mfrc -v MLD,MLD_2mfrc result_2mfrc.nc

ncrename -v MLD,MLD_1mfom result_1mfom.nc
ncrename -v MLD,MLD_2mfom result_2mfom.nc

ncrename -v Qflux,Qflux_1 -v MLD,MLD_1 1.nc
ncrename -v Qflux,Qflux_2 -v MLD,MLD_2 2.nc

ncrename -v Qflux,Qflux_frc -v MLD,MLD_frc frc.nc

ncrename -v MLD,MLD_fom fom.nc


ncks -A result_1m2.nc result.nc
ncks -A result_1mfrc.nc result.nc
ncks -A result_2mfrc.nc result.nc
ncks -A result_1mfom.nc result.nc
ncks -A result_2mfom.nc result.nc

ncks -A 1.nc   -v Qflux_1,MLD_1     result.nc
ncks -A 2.nc   -v Qflux_2,MLD_2     result.nc
ncks -A frc.nc -v Qflux_frc,MLD_frc result.nc
ncks -A fom.nc -v MLD_fom           result.nc

cd ../
mv tmp/result.nc .

rm -rf tmp
