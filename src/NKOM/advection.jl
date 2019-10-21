function calDiffAdv_QUICK!(
    ocn :: Ocean,
)
    calHorVelBnd!(
        Nx    = ocn.Nx,
        Ny    = ocn.Ny,
        Nz    = ocn.Nz,
        weight_e = gi.weight_e,
        weight_n = gi.weight_n,
        u     = ocn.u,
        v     = ocn.v,
        u_bnd = ocn.u_bnd,
        v_bnd = ocn.v_bnd,
        mask3 = ocn.mask3,
    )

    calVerVelBnd!(
        gi = ocn.gi,
        Nx    = ocn.Nx,
        Ny    = ocn.Ny,
        Nz    = ocn.Nz,
        u_bnd = ocn.u_bnd,
        v_bnd = ocn.v_bnd,
        w_bnd = ocn.w_bnd,
        hs    = ocn.hs,
        div   = ocn.div,
        mask3 = ocn.mask3,
    )

    calGRAD!()
    assignHorPhiStar!() 
    assignVerPhiStar!() 
    

end

function calVerVelBnd!(;
    gi       :: GridInfo,
    Nx       :: Integer,
    Ny       :: Integer,
    Nz       :: AbstractArray{Integer, 2},
    u_bnd    :: AbstractArray{Float64, 3},   # ( Nz_bone  , Nx+1, Ny   )
    v_bnd    :: AbstractArray{Float64, 3},   # ( Nz_bone  , Nx  , Ny+1 )
    w_bnd    :: AbstractArray{Float64, 3},   # ( Nz_bone+1, Nx  , Ny   )
    hs       :: AbstractArray{Float64, 3},   # ( Nz_bone  , Nx  , Ny   )
    div      :: AbstractArray{Float64, 3},   # ( Nz_bone  , Nx  , Ny   )
    mask3    :: AbstractArray{Float64, 3},   # ( Nz_bone  , Nx  , Ny   )
)

    for i=1:Nx, j=1:Ny

        w_bnd[1, i, j] = 0.0

        for k=1:Nz[i, j]

            if mask3[k, i, j] == 0.0
                break
            end
            
            flux_w = u_bnd[k, i,   j  ]  * gi.ds4[i, j]
            flux_e = u_bnd[k, i+1, j  ]  * gi.ds2[i, j]

            flux_s = v_bnd[k, i,   j  ]  * gi.ds1[i, j]
            flux_n = v_bnd[k, i,   j+1]  * gi.ds3[i, j]

            div[k, i, j] =  (  
                u_bnd[k, i+1, j  ]  * gi.ds2[i, j]
              - u_bnd[k, i,   j  ]  * gi.ds4[i, j]
              + v_bnd[k, i,   j+1]  * gi.ds3[i, j]
              - v_bnd[k, i,   j  ]  * gi.ds1[i, j]
            ) / gi.dσ[i, j]

            w_bnd[k+1, i, j] = w_bnd[k, i, j] + div[k, i, j] * hs[k, i, j]
        end
    end
 
end



function calHorVelBnd!(;
    Nx       :: Integer,
    Ny       :: Integer,
    Nz       :: AbstractArray{Integer, 2},
    weight_e :: AbstractArray{Float64, 2},   # (Nx+1, Ny)
    weight_n :: AbstractArray{Float64, 2},   # (Nx, Ny+1)
    u        :: AbstractArray{Float64, 3},   # (Nz_bone, Nx, Ny)
    v        :: AbstractArray{Float64, 3},   # (Nz_bone, Nx, Ny)
    u_bnd    :: AbstractArray{Float64, 3},   # (Nz_bone, Nx+1, Ny)
    v_bnd    :: AbstractArray{Float64, 3},   # (Nz_bone, Nx, Ny+1)
    mask3    :: AbstractArray{Float64, 3},   # (Nz_bone, Nx, Ny)
)

    for i=2:Nx, j=1:Ny
        for k=1:Nz[i, j]
            if mask3[k, i, j] == 0.0 || mask3[k, i-1, j] == 0.0
                u_bnd[k, i, j] = 0.0
            else
                u_bnd[k, i, j] = u[k, i-1, j] * (1.0 - weight_e[i, j]) + u[k, i, j] * weight_e[i, j]
            end
        end
    end
    
    for j=1:Ny
        for k=1:Nz[i, j]
            if mask3[k, 1, j] == 0.0 || mask3[k, Nx, j] == 0.0
                u_bnd[k, 1, j] = 0.0
            else
                u_bnd[k, 1, j] = u_bnd[k, Nx+1, j] = u[k, Nx, j] * (1.0 - weight_e[1, j]) + u[k, 1, j] * weight_e[1, j]
            end
        end
    end

    for i=1:Nx, j=2:Ny-1
        for k=1:Nz[i, j]
            if mask3[k, i, j-1] == 0.0 || mask3[k, i, j] == 0.0
                v_bnd[k, i, j] = 0.0
            else
                v_bnd[k, i, j] = v[k, i, j-1] * (1.0 - weight_n[i, j]) + v[k, i, j] * weight_n[i, j]
            end
        end
    end

end

function assignHorPhiStar!(;
)

end

function assignVerPhiStar!(;
)
end

function calHorAdv!(;
)
end

function calVerAdv!(;
)
end
