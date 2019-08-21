mutable struct Configs
    Δt            :: Float64
    substeps      :: Integer
    use_qflx      :: Bool
    use_h_ML      :: Bool
    do_diffusion  :: Bool
    do_relaxation :: Bool
    do_convadjust :: Bool
    rad_scheme    :: Symbol
end
