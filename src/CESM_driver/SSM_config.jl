#=

# Example:
caseroot    = "/home/tienyiah/projects/cesm1_test/CICE_f45"
wdir        = "/home/tienyiah/cesm1/scratch/CICE_f45/run"
domain_file = "/home/tienyiah/cesm_inputdata/cesm1/share/domains/domain.ocn.gx3v7.120323.nc"
zs = collect(Float64, range(0, -500, step=-5))
K = 1e-5
max_try = 60
output_record_length = 365

=#




zs = collect(Float64, range(0, -500, step=-5))
K = 1e-5
max_try = 60
output_record_length = 365
