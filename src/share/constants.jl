

T_fw_frz = 273.15              # Freeze point of freshwater in Kelvin
T_sw_frz = T_fw_frz - 1.8      # Freeze point of seawater in Kelvin

T_ref = T_fw_frz + 20.0        # Reference temperature of thermal expansion coefficient / salinity coefficient
α   = 3e-4     # K^-1    http://www.kayelaby.npl.co.uk/general_physics/2_7/2_7_9.html
β   = 1e-3     # Simple estimation
c_p = 3985.0   # J / kg / K
ρ   = 1027.0   # kg / m^3
g   = 9.8      # m / s^2

αgρc = α * g / (ρ * c_p)

b_sw_frz = α * g * (T_sw_frz - T_ref)

missing_value = 1e20
