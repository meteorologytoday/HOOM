"""

    calWeOrMLD(; h, B, fric_u, Δb, m=0.45, n=0.2)

# Description


 The calculation follows eqn (30) of Gasper (1988). The changes are:

   1. We assume that all sunlight is absorbed at the very surface. According to Paulson and James (1988), in upper ocean
      90% of the radiation is absorbed within the first 10 meters. Since even POP has the resultion of 10m, it is reasonable
      that we ignore the structrue of light absorbtion in sea-water.

   2. Turbulent kinetic energy is added into the equation to avoid the case Δb = 0 which results in we → ∞. This follows
      the reasoning in Alexander et al. (2000)

 References: 

   Gaspar, Philippe. "Modeling the seasonal cycle of the upper ocean." Journal of Physical Oceanography 18.2 (1988): 161-180.

   Paulson, Clayton A., and James J. Simpson. "Irradiance measurements in the upper ocean." 
   Journal of Physical Oceanography 7.6 (1977): 952-956.

   Alexander, Michael A., James D. Scott, and Clara Deser. "Processes that influence sea surface temperature and ocean mixed 
   layer depth variability in a coupled model." Journal of Geophysical Research: Oceans 105.C7 (2000): 16823-16842.


# Return values
This function returns a list with two elements. The first is a symbol. ``:we`` indicates the second value is the entrainment speed whereas ``:MLD`` indicates the second value is the diagnosed MLD.

"""
function calWeOrMLD(;
    h_ML   :: Float64,
    B      :: Float64, 
    fric_u :: Float64,  
    Δb     :: Float64,
    f      :: Float64,
    m::Float64 = 0.8,
    n::Float64 = 0.20,
)


    if Δb < 0
        throw(ErrorException("Δb cannot be negative. Right now Δb = ", Δb))
    end

    Term1 = 2.0 * m * fric_u^3.0 * exp( -h_ML / (fric_u / abs(f)))
    Term2 = 0.5 * (B * (1.0 + n) - abs(B) * (1.0 - n))

    RHS = Term1 + h_ML * Term2

    if RHS > 0
        k = getTKE(fric_u=fric_u)
        we = RHS / (h_ML * Δb + k)
        #if we > 1e-3
        #    println("we abnormally large: ", we)
        #end
        #println(":we, h: ", h, "; Δb: ", Δb, "; B: ", B, "; k:", k)
        return :we, we
    else

        # h becomes diagnostic.
       
        if Term2 == 0
            h_ML_diag = h_ML
        else
            h_ML_diag = - Term1 / Term2
        end
    
        return :MLD, h_ML_diag
    end

end


