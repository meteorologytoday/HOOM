using Pkg 


plist = "Formatting, NCDatasets, StatsBase, ArgParse"

for pname in split(replace(plist, " "=>""), ",")
    println("Target Package: ", pname)
    Pkg.add(String(pname))
end

