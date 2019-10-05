To prescribe sea-ice,

- Use CICE%PRES mode.
- Setup "SSTICE_DATA_FILENAME", "SSTICE_GRID_FILENAME", "SSTICE_YEAR_ALIGN", "SSTICE_YEAR_START", "SSTICE_YEAR_END" in env_run.xml. Year values are 1 based. That is, the first year is year 1.
- Reading from `cesm1221/models/ice/cice/src/drivers/cpl_share/ice_prescribed_mod.F90`, variable names of forcing file : (1) time (the only 1D variables) (2) ice_conv : seaice-fraction, (3) lat : latitude, (4) lon : longitude, (5) area (6) mask
- Notice that seaice thickness is being prescribed as 2m in NH and 1m in SH. Also seaice thickness is set to 0m within 40S-40N (`ice_prescribed_mod.F90` line 455-463).

