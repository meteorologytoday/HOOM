function writeLog(args...; force :: Bool = true)
    
    comm = MPI.COMM_WORLD
    rank = MPI.Comm_rank(comm)

    if force || rank == 0
        println(format(args...))
    end
end
