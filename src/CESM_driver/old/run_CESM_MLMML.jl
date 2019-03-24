using ArgParse
using Printf

include("../default_config.jl")

# ===== Load Core Module =====
include("CESM_CORE_MLMML.jl")
OMMODULE = Main.CESM_CORE_MLMML

# ===== Load General Driver =====
include("../driver.jl")
