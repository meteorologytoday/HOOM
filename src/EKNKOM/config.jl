mutable struct Config
    substeps      :: Integer
    use_qflx      :: Bool
    use_h_ML      :: Bool
    Δt            :: Float64
    do_diffusion  :: Bool
    do_relaxation :: Bool
    do_convadjust :: Bool
    rad_scheme    :: Symbol
end
