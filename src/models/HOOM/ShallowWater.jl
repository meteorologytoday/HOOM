
module ShallowWaterModule

    mutable struct SWStatus
        
        Nx :: Integer
        Ny :: Integer
        Nz :: Integer

        z_bnd :: AbstractArray{Float64, 1}
        
        u  :: AbstractArray{Float64, 3}
        v  :: AbstractArray{Float64, 3}
        T  :: AbstractArray{Float64, 3} 
        S  :: AbstractArray{Float64, 3} 
        X  :: AbstractArray{Float64, 3}   # passive tracer
 
        # barotropic velocity
        U  :: AbstractArray{Float64, 2}
        V  :: AbstractArray{Float64, 2}
        η  :: AbstractArray{Float64, 2}

       
        function SWStatus(
            Nx :: Integer,
            Ny :: Integer,
            Nz :: Integer,
            z_bnd :: AbstractArray{Float64, 1}
        )
       

            u = zeros(Float64, Nx, Ny, Nz)
            v = zeros(Float64, Nx, Ny, Nz)
            b = zeros(Float64, Nx, Ny, Nz)

            U = zeros(Float64, Nx, Ny)
            V = zeros(Float64, Nx, Ny)
            η = zeros(Float64, Nx, Ny)


            return new(
                Nx, Ny, Nz,
                z_bnd,
                u, v, b, U, V, η
            ) 
        end

    end

    mutable struct ABIIIObj
        gi :: GridInfo
        
    
        function ABIIIObj(
            ShallowWater
        )
                        
        end
    end


    function ABIII!(
        o :: ABIIIObj
    )

    end
end
