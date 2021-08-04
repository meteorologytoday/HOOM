mutable struct Env

    config :: Dict
    
    sub_yrng    :: Any
    Nx :: Integer
    Ny :: Integer
    Nz :: Integer

    function Env(
        config;
        sub_yrng = nothing,
    )
        
        validated_config = validateConfig(config, genConfigEntryList(:ENV)) 
           
        # mask =>   lnd = 0, ocn = 1
        gf = PolelikeCoordinate.CurvilinearSphericalGridFile(
            validated_config[:domain_file];
            R  = 6371229.0,
            Ω  = 2π / (86400 / (1 + 365/365)),
        )

        Nx = gf.Nx
        Ny = (sub_yrng == nothing) ? gf.Ny : length(sub_yrng)
        Nz = length(validated_config[:z_w]) - 1

        return new(
            
            validated_config,
 
            sub_yrng,        
            Nx,
            Ny,
            Nz,

        )
    end

end


