merge!(overwrite_configs, Dict(
    :MLD_scheme                   => :datastream,
    :Qflux_scheme                 => :energy_flux,
    :vertical_diffusion_scheme    => :off,
    :horizontal_diffusion_scheme  => :off,
    :relaxation_scheme            => :on,
    :convective_adjustment_scheme => :on,
    :radiation_scheme             => :step,
    :advection_scheme             => :ekman_codron2012_partition,
))
