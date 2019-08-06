using JSON
using Formatting

target_year = 200

archive_dir = "~/scratch-tienyiao/archive/"
cases_dir = "./cases"
all_casenames = split(read(pipeline(`ls $cases_dir`), String), "\n", keepempty=false)
working_casenames = split(read(pipeline(`qstat -f`, `grep Job_Name`, `awk '{print $3}'`), String), "\n", keepempty=false)

println("##### All casenames: #####")
JSON.print(all_casenames, 4)

println("##### Working casenames: #####")
JSON.print(working_casenames, 4)


failed_casenames = []
for (i, casename) in enumerate(all_casenames)

    if casename in working_casenames
        continue
    end

    target_file = joinpath(
        archive_dir,
        casename,
        "ocn", "hist",
        format("{:s}.ocn.h.monthly.{:04d}.nc", casename, target_year)
    )

    if !isfile(target_file)
        push!(failed_casenames, casename)
    end
end


println("##### Failed casenames: #####")
JSON.print(failed_casenames, 4)

working_dir = pwd()
println("Do you want to qsub these cases? (yes/no)")
if strip(readline(stdin)) == "yes"
    for casename in failed_casenames
        cd(working_dir)
        println("Submitting ", casename)
        cd(joinpath(cases_dir, casename))
        run(`qsub $(casename).run`)
    end
end

println("##### END #####")
