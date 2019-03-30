using Distributed


addprocs(2)
println(procs())


mod_name = "moduleA"
mod_file = mod_name * ".jl"
mod_sb = Symbol(mod_name)
for p in procs()
    println("Doing worker: ", p)
    remotecall_fetch(include, p, mod_file)
    remotecall_fetch(eval, p, :(using .$mod_sb))
    remotecall_fetch(eval, p, :($mod_sb.funA())) 
end



