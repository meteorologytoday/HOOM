using Pkg 


plist = "Formatting, NCDatasets, StatsBase, ArgParse"

println("Going to install the following pkgs: " * plist)
for pname in split(replace(plist, " "=>""), ",")
    println("Target Package: ", pname)
    Pkg.add(String(pname))
end

