function genColView(
    arr :: Union{AbstractArray{T, 3}, AbstractArray{T, 4}}
) where T

    # vectorization the for loop is not faster

    s = size(arr)

    if length(s) == 3

        _, Nx, Ny = size(arr)
        view_arr  = Array{SubArray}(undef, Nx, Ny)

        for i=1:Nx, j=1:Ny
            view_arr[i, j] = view(arr, :, i, j)
        end

    elseif length(s) == 4

        _, Nx, Ny, NX = size(arr)
        view_arr  = Array{SubArray}(undef, Nx, Ny, NX)

        for i=1:Nx, j=1:Ny, x =1:NX
            view_arr[i, j, x] = view(arr, :, i, j, x)
        end

    end

    return view_arr
end

