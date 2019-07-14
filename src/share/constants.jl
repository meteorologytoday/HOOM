const T_fw_frz = 273.15              # Freeze point of freshwater in Kelvin
const T_sw_frz = T_fw_frz - 1.8      # Freeze point of seawater in Kelvin

const T_ref = T_fw_frz + 20.0        # Reference temperature of thermal expansion coefficient / salinity coefficient
const S_ref = 35.0                   # Reference salinity    of thermal expansion coefficient / salinity coefficient  (PSU)
const c_p = 3996.0   # J / kg / K   copied from models/csm_share/shr/shr_const_mod.F90
const ρ   = 1026.0   # kg / m^3     copied from models/csm_share/shr/shr_const_mod.F90
const g   = 9.80616  # m / s^2      copied from models/csm_share/shr/shr_const_mod.F90
const α   = 3e-4     # K^-1    http://www.kayelaby.npl.co.uk/general_physics/2_7/2_7_9.html
const β   = 1e-3     # Simple estimation

const αgρc = α * g / (ρ * c_p)

const b_sw_frz = α * (T_sw_frz - T_ref)


const S_surf_avg = 35.0    # The average surface salinity used to update salinity when raining/evaporating.


const missing_value = 1e20

const Re = 6371e3          # Radius of Earth
const Ωe = 2.0 * π * (1.0 + 1.0/365.0) / 86400.0
