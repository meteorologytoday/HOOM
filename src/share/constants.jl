const T_fw_frz_Kelvin = 273.15

const T_fw_frz =  0.0       # Freeze point of freshwater in Celcius
const T_sw_frz = -1.8       # Freeze point of seawater in Celcius

const c_p   = 3996.0   # J / kg / C   copied from models/csm_share/shr/shr_const_mod.F90
const ρ     = 1026.0   # kg / m^3     copied from models/csm_share/shr/shr_const_mod.F90
const ρ_fw  = 1000.0   # kg / m^3     copied from models/csm_share/shr/shr_const_mod.F90
const g     = 9.80616  # m / s^2      copied from models/csm_share/shr/shr_const_mod.F90

const ρc = ρ * c_p

const missing_value = 1e20

const Re = 6371e3          # Radius of Earth
const Ωe = 2.0 * π * (1.0 + 1.0/365.0) / 86400.0
