function OC_doConvectiveAdjustment!(
        occ :: OceanColumnCollection,
        i   :: Integer,
        j   :: Integer,
    )

    if_adjust, occ.T_ML[i, j], occ.S_ML[i, j], occ.h_ML[i, j], occ.FLDO[i, j] = doConvectiveAdjustment!(
        zs       = occ.zs_vw[i, j],
        bs       = occ.bs_vw[i, j],
        Ts       = occ.Ts_vw[i, j],
        Ss       = occ.Ss_vw[i, j],
        h_ML     = occ.h_ML[i, j],
        b_ML     = occ.b_ML[i, j],
        T_ML     = occ.T_ML[i, j],
        S_ML     = occ.S_ML[i, j],
        FLDO     = occ.FLDO[i, j],
        Nz       = occ.Nz[i, j],
        h_ML_max = occ.h_ML_max[i, j],
    )

    return if_adjust
end


"""

This function only do convective adjustment for the upper most mixed layer.
It searches for the lowest layer that has larger buoyancy than mixed-layer then mixed all layers above it.

By default it only mixes T and S but not b


"""
function doConvectiveAdjustment!(;
    zs   :: AbstractArray{Float64, 1},
    bs   :: AbstractArray{Float64, 1},
    Ts   :: AbstractArray{Float64, 1},
    Ss   :: AbstractArray{Float64, 1},
    h_ML :: Float64,
    b_ML :: Float64,
    T_ML :: Float64,
    S_ML :: Float64,
    FLDO :: Integer,
    Nz   :: Integer,
    h_ML_max :: Float64,
)
    
    if_adjust = false

    if FLDO == -1
        return if_adjust, T_ML, S_ML, h_ML, FLDO 
    end

    # 1. Search from bottom to see if buoyancy is monotically increasing
    # 2. If not, record the peak, then keep detecting until hitting the 
    #    layer X_top having equal or greater buoyancy. Record this interval.
    # 3. Find the minimum value in this interval b_min.
    # 4. Use b_min to decide the bottom layer X_bot going to be mixed (Search 
    #    downward).
    # 5. Mix layers between X_bot ~ X_top.

    new_T_ML = T_ML
    new_S_ML = S_ML
    new_h_ML = h_ML
    new_FLDO = FLDO

    stage = :reset
    peak_layer = 0
    top_layer = 0
    bot_layer = 0
    b_peak = 0.0



   for i = Nz:-1:FLDO


        if stage == :reset
            peak_layer = 0
            top_layer = 0
            bot_layer = 0
            b_peak = 0.0
            stage = :search_peak_layer
        end

        if stage == :search_peak_layer

            Δb = bs[i] - ((i==FLDO) ? b_ML : bs[i-1])
            #println("FLDO:", FLDO, "; i:", i, "; Δb:", Δb)
            if Δb > 0.0
                if_adjust = true
                stage = :search_top_layer
                peak_layer = i
                b_peak = bs[peak_layer]
            else
                continue
            end
        end

        if stage == :search_top_layer

            #println(":search_top_layer")
            if i == FLDO
                top_layer = (b_ML > b_peak) ? FLDO : -1
                stage = :search_bot_layer
            elseif bs[i-1] > b_peak
                top_layer = i
                stage = :search_bot_layer
            else
                continue
            end
        end

        if stage == :search_bot_layer

            #println(":search_bot_layer")

            if peak_layer == Nz

                bot_layer = peak_layer
                stage = :start_adjustment

            else
                b_min = 0.0
                if top_layer == -1
                    b_min = min(b_ML, minimum(bs[FLDO:peak_layer]))
                else
                    b_min = minimum(bs[top_layer:peak_layer])
                end 

                bot_layer = peak_layer + 1
                while true
                    if bs[bot_layer] >= b_min
                        if bot_layer == Nz
                            stage = :start_adjustment
                            break
                        else
                            bot_layer += 1
                        end
                    else
                        stage = :start_adjustment
                        break
                    end
                end
            end 
        end


        if stage == :start_adjustment
            #println(":start_adjustment")

            
            bot_z = zs[bot_layer+1]
            #println(zs[bot_layer+1]," v.s.  ", -h_ML_max, "; bot_z: ", bot_z)
            top_z = (top_layer == -1) ? 0.0 : (
                 (top_layer == FLDO) ? -h_ML : zs[top_layer]
            )
            Δz = top_z - bot_z

            mixed_T = (getIntegratedQuantity(
                zs       =  zs,
                qs       =  Ts,
                q_ML     =  T_ML,
                h_ML     =  h_ML,
                Nz       =  Nz,
                target_z =  bot_z
            ) - getIntegratedQuantity(
                zs       =  zs,
                qs       =  Ts,
                q_ML     =  T_ML,
                h_ML     =  h_ML,
                Nz       =  Nz,
                target_z =  top_z
            ))  / Δz
 
            mixed_S = (getIntegratedQuantity(
                zs       =  zs,
                qs       =  Ss,
                q_ML     =  S_ML,
                h_ML     =  h_ML,
                Nz       =  Nz,
                target_z =  bot_z
            ) - getIntegratedQuantity(
                zs       =  zs,
                qs       =  Ss,
                q_ML     =  S_ML,
                h_ML     =  h_ML,
                Nz       =  Nz,
                target_z =  top_z
            ))  / Δz
           
            if top_layer == -1  # Even the mixed layer is mixed

                # 2019/08/04 Decide that convective adjustment does not change MLD.
                # This makes ML dynamic less complicated

                new_T_ML = mixed_T 
                new_S_ML = mixed_S
                
                # update T, S profile but do not update h_ML and FLDO 
                setMixedLayer!(Ts=Ts, Ss=Ss, zs=zs, T_ML=new_T_ML, S_ML=new_S_ML, h_ML= - bot_z, Nz=Nz)

                #= 

                new_T_ML = mixed_T
                new_S_ML = mixed_S
                new_h_ML = - bot_z

                new_FLDO = setMixedLayer!(Ts=Ts, Ss=Ss, zs=zs, T_ML=new_T_ML, S_ML=new_S_ML, h_ML=new_h_ML, Nz=Nz)
                
                if new_h_ML > h_ML_max

                    # modified on 2019/05/10
                    new_h_ML = h_ML_max
                    new_FLDO = setMixedLayer!(Ts=Ts, Ss=Ss, zs=zs, T_ML=new_T_ML, S_ML=new_S_ML, h_ML=new_h_ML, Nz=Nz)

                    # original code before 2019/05/10
                    # new_FLDO = getFLDO(zs=zs, h_ML=h_ML_max, Nz=Nz)   # original code
                end
                =#

            else
                Ts[top_layer:bot_layer] .= mixed_T
                Ss[top_layer:bot_layer] .= mixed_S
            end 

            stage = :reset 
        end

    end

#    if if_adjust
#        println("ADJUST! ", h_ML, " => ", new_h_ML)
#    end

    return if_adjust, new_T_ML, new_S_ML, new_h_ML, new_FLDO
end

           

