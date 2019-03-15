"""

    stepOceanColumn!(;
        oc      :: OceanColumn,
        u_fric  :: Float64,
        B0      :: Float64,
        J0      :: Float64,
        Δt      :: Float64 
    )

# Description
This function update the OceanColumn forward in time.

"""
function stepOceanColumnCollection!(
    occ    :: OceanColumnCollection;
    fric_u :: AbstractArray{Float64, 2}, # Currently assumed to be u10
    B0     :: AbstractArray{Float64, 2},
    J0     :: AbstractArray{Float64, 2},
    Δt     :: Float64 
)

    for i = 1:occ.Nx, j = 1:occ.Ny

        if occ.mask[i, j] == 0.0
            continue
        end

        # Pseudo code
        # Current using only Euler forward scheme:
        # 1. Determine h at t+Δt
        # 2. Determine how many layers are going to be
        #    taken away by ML.
        # 3. Cal b at t+Δt for both ML and DO
        # 4. Detect if it is buoyantly stable.
        #    Correct it (i.e. convection) if it is not.
        # 5. If convection happens, redetermine h.

        # p.s.: Need to examine carefully about the
        #       conservation of buoyancy in water column

        #println("### h: ", oc.h)
        #println("FLDO:", oc.FLDO)

        FLDO = occ.FLDO[i, j]
        Δb = (FLDO != -1) ? occ.b_ML[i, j] - occ.bs[i, j, FLDO] : 0.0

        # After convective adjustment, there still might
        # be some numerical error making Δb slightly negative
        # (the one I got is like -1e-15). So I set a tolarence
        # ( 0.001 K ≈ 3e-6 m/s^2 ).
        if Δb < 0.0 && abs(Δb) <= 3e-6
            Δb = 0.0
        end

        total_flx = B0[i, j] + J0[i, j]

        #fric_u = getFricU(ua=ua)
        flag, val = calWeOrMLD(;
            h_ML   = occ.h_ML[i, j],
            B      = total_flx,
            fric_u = fric_u[i, j],
            Δb     = Δb
        )
        #println("Before:" , oc.bs[10], "; oc.FLDO = ", oc.FLDO, "; Δb = ", Δb)
        #println("B: ", total_flx , ";Δb: ", Δb , "; fric_u: ", fric_u[i, j])

        # 1
        if flag == :MLD
            we = 0.0
            new_h_ML  = val
        elseif flag == :we
            we = val 
            new_h_ML = occ.h_ML[i, j] + Δt * we
        end
        new_h_ML = boundMLD(new_h_ML; h_ML_max=min(h_ML_max, occ.zs[1] - occ.zs[end]))

        #println("flag: ", String(flag), "; val: ", val, "; new_h_ML: ", new_h_ML)
        # 2
        # 3

        # ML
        #      i: Calculate integrated buoyancy that should
        #         be conserved purely through entrainment
        #     ii: Add to total buoyancy

        hb_new = OC_getIntegratedBuoyancy(occ, i, j; target_z = -new_h_ML)
      
        hb_chg_by_F = - total_flx * Δt

        #println(new_h, "; ", hb_new, ", ")
        new_b_ML = (hb_new + hb_chg_by_F) / new_h_ML
        
        OC_setMixedLayer!(
            occ, i, j;
            b_ML=new_b_ML,
            h_ML=new_h_ML,
        )
        
        OC_doDiffusion_EulerBackward!(occ, i, j; Δt=Δt)
        OC_doConvectiveAdjustment!(occ, i, j)

        # Freeze potential. Calculation mimics the one written in CESM1 docn_comp_mod.F90
        occ.qflx2atm[i, j] = (T_sw_frz - b2T(occ.b_ML[i, j])) * ρ * c_p * occ.h_ML[i, j] / Δt
        occ.b_ML[i, j] = max(b_sw_frz, occ.b_ML[i, j])
        
    end

    updateSST!(occ)
end


