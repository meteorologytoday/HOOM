const T_fw_frz = 273.15              # Freeze point of freshwater in Kelvin
const T_sw_frz = T_fw_frz - 1.8      # Freeze point of seawater in Kelvin

const T_ref = T_fw_frz + 20.0        # Reference temperature of thermal expansion coefficient / salinity coefficient
const S_ref = 35.0                   # Reference salinity    of thermal expansion coefficient / salinity coefficient  (PSU)
const c_p = 3996.0   # J / kg / K   copied from models/csm_share/shr/shr_const_mod.F90
const ρ   = 1026.0   # kg / m^3     copied from models/csm_share/shr/shr_const_mod.F90
const g   = 9.80616  # m / s^2      copied from models/csm_share/shr/shr_const_mod.F90


#
# The numbers are referenced from Bryan and Cox (1972) at depth Z = 250m. 
#
# Formulas for α and β are:
#
#    α := - (1/ρ)(∂ρ/∂T) 
#    β :=   (1/ρ)(∂ρ/∂S) 
#
# Reference:
#
#     Bryan, Kirk, and Michael D. Cox. "An approximate equation of state for numerical models of ocean circulation." Journal of Physical Oceanography 2.4 (1972): 510-514.
#
const α   = 1.5781e-4 / 1.028475     # 1/K
const β   = 7.8318e-4 / 1.028475     # 1/PSU

const αgρc = α * g / (ρ * c_p)

const b_sw_frz = α * (T_sw_frz - T_ref)
const S_surf_avg = 35.0    # The average surface salinity used to update salinity when raining/evaporating.


const missing_value = 1e20

const Re = 6371e3          # Radius of Earth
const Ωe = 2.0 * π * (1.0 + 1.0/365.0) / 86400.0
