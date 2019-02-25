using Pkg 


plist = "Formatting, NCDatasets, StatsBase, Hwloc"

for pname in split(replace(plist, " "=>""), ",")
    println("Target Package: ", pname)
    Pkg.add(String(pname))
end

