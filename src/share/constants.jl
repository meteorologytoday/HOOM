const T_fw_frz_Kelvin = 273.15

const T_fw_frz =  0.0       # Freeze point of freshwater in Celcius
const T_sw_frz = -1.8       # Freeze point of seawater in Celcius

const T_ref = 20.0          # Reference temperature of thermal expansion coefficient / salinity coefficient in Celcius
const S_ref = 35.0          # Reference salinity    of thermal expansion coefficient / salinity coefficient  (PSU)
const c_p   = 3996.0   # J / kg / C   copied from models/csm_share/shr/shr_const_mod.F90
const ρ     = 1026.0   # kg / m^3     copied from models/csm_share/shr/shr_const_mod.F90
const ρ_fw  = 1000.0   # kg / m^3     copied from models/csm_share/shr/shr_const_mod.F90
const g     = 9.80616  # m / s^2      copied from models/csm_share/shr/shr_const_mod.F90

# Follow Bryan and Cox (1972) : An Approximate Equation of State for Numerical Models of Ocean Circulation
# Table 2: Take mean of Z = 0,250,500 m. Note that the coefficient X1 and X2 need to be divided by 1000
#          and equation (1) is wrong because 1000 should be dividing. Also note that the signed is reversed
#          because we are talking about buoyancy instead of density.

const α     = 1.633433e-4     # per degC
const β     = 7.814766e-4     # per PSU

const αgρc = α * g / (ρ * c_p)
const ρc = ρ * c_p

const b_sw_frz = α * (T_sw_frz - T_ref)


const S_surf_avg = 35.0    # The average surface salinity used to update salinity when raining/evaporating.


const missing_value = 1e20

const Re = 6371e3          # Radius of Earth
const Ωe = 2.0 * π * (1.0 + 1.0/365.0) / 86400.0
