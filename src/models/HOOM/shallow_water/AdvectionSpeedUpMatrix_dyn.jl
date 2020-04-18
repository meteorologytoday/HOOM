
using SparseArrays

@inline function flat_i(
    k :: Int64,
    i :: Int64,
    j :: Int64,
    Nz :: Int64,
    Nx :: Int64,
    Ny :: Int64,
)
    return k + (i-1) * Nz + (j-1) * Nz * Nx
end

@inline function cyc(i::Int64, N::Int64)
    return mod(i-1, N) + 1
end

@inline function speye(dtype, n)
    return spdiagm(0=>ones(dtype, n))
end

@inline function dropall!(a)
    replace!(a, NaN=>0)
    dropzeros!(a)
end


# Assuming x-direction is periodic
struct MatrixOperators

    # Nomenclature:
    #
    # [new-grid][direction][old-grid]
    #
    # U_W_T : sending variable westward from T grid to U grid

    U_I_U
    V_I_V
    T_I_T

    U_W_T
    U_E_T
    U_W_U
    U_E_U
    V_W_V
    V_E_V

    V_S_T
    V_N_T
    U_S_U
    U_N_U
    V_S_V
    V_N_V

    T_S_V
    T_N_V
    T_W_U
    T_E_U 
    
    U_SW_V
    U_SE_V
    U_NW_V
    U_NE_V
    V_SW_U
    V_SE_U
    V_NW_U
    V_NE_U


    #= z related
    T_UP_T
    T_DN_T
    T_UP_W
    T_DN_W
    W_UP_T
    W_DN_T
    =#

    function MatrixOperators(;
        Nx             :: Int64,
        Ny             :: Int64,
        Nz             :: Int64,
    )

        # Making operator
        U_dim = (Nz, Nx, Ny)
        V_dim = (Nz, Nx, Ny+1)
        W_dim = (Nz+1, Nx, Ny)
        T_dim = (Nz, Nx, Ny)

        U_pts = reduce(*, U_dim)
        V_pts = reduce(*, V_dim)
        W_pts = reduce(*, W_dim)
        T_pts = reduce(*, T_dim)

        U_I_U = speye(Float64, U_pts)
        V_I_V = speye(Float64, V_pts)
        W_I_W = speye(Float64, W_pts)
        T_I_T = speye(Float64, T_pts)

        U_I_U_expand = vcat(U_I_U, zeros(Float64, 1, U_pts))
        V_I_V_expand = vcat(V_I_V, zeros(Float64, 1, V_pts))
        W_I_W_expand = vcat(W_I_W, zeros(Float64, 1, W_pts))
        T_I_T_expand = vcat(T_I_T, zeros(Float64, 1, T_pts))


        num_U = zeros(Int64, U_dim...)
        num_V = zeros(Int64, V_dim...)
        num_W = zeros(Int64, W_dim...)
        num_T = zeros(Int64, T_dim...)

        num_U[:] = 1:length(num_U)
        num_V[:] = 1:length(num_V)
        num_W[:] = 1:length(num_W)
        num_T[:] = 1:length(num_T)
        
        U = num_U * 0
        V = num_V * 0
        W = num_W * 0
        T = num_T * 0

        function build!(id_mtx, idx; wipe=:none)
            local result
            rows = size(id_mtx)[1]
            if wipe == :n
                idx[:, :, end] .= rows
            elseif wipe == :s
                idx[:, :,   1] .= rows
            elseif wipe != :none
                throw(ErrorException("Wrong keyword"))
            end

            result = id_mtx[view(idx, :), :]
            dropzeros!(result)
            idx .= 0 # clean so that debug is easir when some girds are not assigned
            return result
        end

        # east, west passing mtx
        U[:, :, :] = num_T;                             U_W_T = build!(T_I_T, U)
        U[:, :, :] = circshift(num_T, (0, 1, 0));       U_E_T = build!(T_I_T, U)

        U[:] = circshift(num_U, (0, -1, 0));            U_W_U = build!(U_I_U, U)
        U[:] = circshift(num_U, (0,  1, 0));            U_E_U = build!(U_I_U, U)

        V[:] = circshift(num_V, (0, -1, 0));            V_W_V = build!(V_I_V, V)
        V[:] = circshift(num_V, (0,  1, 0));            V_E_V = build!(V_I_V, V)
 
        # north, south passing mtx
        V[:, :, 1:Ny]   = num_T;                        V_S_T = build!(T_I_T_expand, V; wipe=:n)
        V[:, :, 2:Ny+1] = num_T;                        V_N_T = build!(T_I_T_expand, V; wipe=:s)

        U[:, :, 1:Ny-1] = view(num_U, :, :, 2:Ny);      U_S_U = build!(U_I_U_expand, U; wipe=:n)
        U[:, :, 2:Ny] = view(num_U, :, :, 1:Ny-1);      U_N_U = build!(U_I_U_expand, U; wipe=:s)
        
        V[:, :, 1:Ny]   = view(num_V, :, :, 2:Ny+1);    V_S_V = build!(V_I_V_expand, V; wipe=:n)
        V[:, :, 2:Ny+1] = view(num_V, :, :, 1:Ny  );    V_N_V = build!(V_I_V_expand, V; wipe=:s)


        # inverse directions
        T_S_V = V_N_T' |> sparse
        T_N_V = V_S_T' |> sparse
        T_W_U = U_E_T' |> sparse
        T_E_U = U_W_T' |> sparse

        # diagonal passing mtx        
        U_SW_V = U_W_T * T_S_V
        U_SE_V = U_E_T * T_S_V
        U_NW_V = U_W_T * T_N_V
        U_NE_V = U_E_T * T_N_V

        V_SW_U = V_S_T * T_W_U
        V_SE_U = V_S_T * T_E_U
        V_NW_U = V_N_T * T_W_U
        V_NE_U = V_N_T * T_E_U

        return new(

            U_I_U, V_I_V, T_I_T,

            U_W_T, U_E_T,
            U_W_U, U_E_U,
            V_W_V, V_E_V,

            V_S_T, V_N_T,
            U_S_U, U_N_U,
            V_S_V, V_N_V,

            T_S_V, T_N_V,
            T_W_U, T_E_U,

            U_SW_V, U_SE_V,
            U_NW_V, U_NE_V,
            V_SW_U, V_SE_U,
            V_NW_U, V_NE_U,

            #T_UP_T, T_DN_T,
            #T_UP_W, T_DN_W,
            #W_UP_T, W_DN_T,
         
        )

    end
end

mutable struct DynamicAdvSpeedUpMatrix
    
    op       :: MatrixOperators

    # Nomenclature:
    #
    # [new-grid][function][old-grid]
    #
    # U_W_T : sending variable westward from T grid to U grid

    U_∂x_T   :: AbstractArray{Float64, 2}   # ∇b 
    V_∂y_T   :: AbstractArray{Float64, 2}   
 
    U_∂x_U   :: AbstractArray{Float64, 2}   # ∇u
    U_∂y_U   :: AbstractArray{Float64, 2}

    V_∂x_V   :: AbstractArray{Float64, 2}   # ∇v
    V_∂y_V   :: AbstractArray{Float64, 2}

    T_DIVx_U :: AbstractArray{Float64, 2}
    T_DIVy_V :: AbstractArray{Float64, 2}
 
    T_Lap_T  :: AbstractArray{Float64, 2}   # Need to solve Poisson equation

    U_interp_V :: AbstractArray{Float64, 2}  # interpolation of V grid onto U grid
    V_interp_U :: AbstractArray{Float64, 2}  # interpolation of U grid onto V grid

    U_f_V      :: AbstractArray{Float64, 2}   # used to get fv on U grid
    V_f_U      :: AbstractArray{Float64, 2}   # used to get fu on V grid
    
    filter_T      :: AbstractArray{Float64, 2}
    filter_U      :: AbstractArray{Float64, 2}
    filter_V      :: AbstractArray{Float64, 2}

    function DynamicAdvSpeedUpMatrix(;
        gi             :: PolelikeCoordinate.GridInfo,
        Nz             :: Int64,
        mask2          :: AbstractArray{Float64, 2},
    )
        Nx = gi.Nx
        Ny = gi.Ny
        op = MatrixOperators(Nx=Nx, Ny=Ny, Nz=Nz)

        # making filter
        # define a converter to make 2D variable repeat in z direction for Nz times
        cvt23 = (x,) -> repeat(reshape(x, 1, size(x)...), outer=(Nz, 1, 1))
        cvt23_diagm = (x,) -> spdiagm( 0 => view(cvt23(x), :) )

        mask3_flat = view(mask2 |> cvt23, :)
        
        if_mask_north_V = op.V_S_T * mask3_flat
        if_mask_south_V = op.V_N_T * mask3_flat
        if_mask_east_U  = op.U_W_T * mask3_flat
        if_mask_west_U  = op.U_E_T * mask3_flat


        #println(size(op.V_S_T))
        #println(size(mask3_flat))

        filter_T = spdiagm(0 => mask3_flat)
        filter_V = spdiagm(0 => if_mask_north_V .* if_mask_south_V)
        filter_U = spdiagm(0 => if_mask_east_U  .* if_mask_west_U)

        #println("if_mask_north_V", size(if_mask_north_V))
        #println("filter_V", size(filter_V))

        dropzeros!(filter_V)
        dropzeros!(filter_U)

        # Some face area and lengths

        U_invΔx_U = (gi.dx_w.^(-1) |> cvt23_diagm)
        U_invΔy_U = (gi.DY.^(-1)   |> cvt23_diagm)
        
        V_invΔx_V = (gi.DX.^(-1)   |> cvt23_diagm)
        V_invΔy_V = (hcat(gi.dy_s, gi.dy_n[:, end]).^(-1) |> cvt23_diagm)

        V_Δx_V    = (gi.DX       |> cvt23_diagm)
        U_Δy_U    = (gi.DY       |> cvt23_diagm)
        T_invΔσ_T = (gi.dσ.^(-1) |> cvt23_diagm)

        # MAGIC!!
        U_∂x_T = filter_U * U_invΔx_U * (op.U_W_T - op.U_E_T) ; dropzeros!(U_∂x_T);

        #println(size(filter_V))
        #println(size(V_invΔy_V))
        #println(size(op.V_S_T))
        #println(size(op.V_N_T))


        V_∂y_T = filter_V * V_invΔy_V * (op.V_S_T - op.V_N_T) ; dropzeros!(V_∂y_T);

        U_∂x_U = filter_U * U_invΔx_U * (op.U_W_U - op.U_E_U) / 2.0 ; dropzeros!(U_∂x_U);
        U_∂y_U = filter_U * U_invΔy_U * (op.U_S_U - op.U_N_U) / 2.0 ; dropzeros!(U_∂y_U);

        V_∂x_V = filter_V * V_invΔx_V * (op.V_W_V - op.V_E_V) / 2.0 ; dropzeros!(V_∂x_V);
        V_∂y_V = filter_V * V_invΔy_V * (op.V_S_V - op.V_N_V) / 2.0 ; dropzeros!(V_∂y_V);
 
        T_DIVx_U = filter_T * T_invΔσ_T * ( op.T_W_U - op.T_E_U  ) * U_Δy_U  ; dropzeros!(T_DIVx_U);
        T_DIVy_V = filter_T * T_invΔσ_T * ( op.T_S_V - op.T_N_V  ) * V_Δx_V  ; dropzeros!(T_DIVy_V);

        U_interp_V = (op.U_SW_V + op.U_SE_V + op.U_NW_V + op.U_NE_V) * filter_V     
        U_interp_V ./= U_interp_V  # weighting. You need to think about it
        dropall!(U_interp_V)

        V_interp_U = (op.V_SW_U + op.V_SE_U + op.V_NW_U + op.V_NE_U) * filter_U
        V_interp_U ./= V_interp_U
        dropall!(V_interp_U)
       
        V_interp_T = (op.V_S_T + op.V_N_T) * filter_T
        V_interp_T ./= V_interp_T
        dropall!(V_interp_T)

        U_interp_T = (op.U_W_T + op.U_E_T) * filter_T
        U_interp_T ./= U_interp_T
        dropall!(U_interp_T)
 
        f = 2.0 * gi.Ω * sin.(gi.c_lat)
        U_f_U = filter_U * spdiagm( 0 => view(U_interp_T * view(f, :), :) )
        V_f_V = filter_V * spdiagm( 0 => view(V_interp_T * view(f, :), :) )
        
        # imagine term fv act on U grid
        # filter_U and filter_V have already been applied in U_f_U and V_f_V
        U_f_V = U_f_U * U_interp_V
        V_f_U = V_f_V * V_interp_U
        

        # Laplacian
        T_Lap_T = T_DIVx_U * U_∂x_T + T_DIVy_V * V_∂y_T 

        return new(
            op,
            U_∂x_T, V_∂y_T,
            U_∂x_U, U_∂y_U,
            V_∂x_V, V_∂y_V,
            T_DIVx_U, T_DIVy_V,
            T_Lap_T,
            U_interp_V, V_interp_U,
            U_f_V, V_f_U,
            filter_T, filter_U, filter_V
        )
    end
end
