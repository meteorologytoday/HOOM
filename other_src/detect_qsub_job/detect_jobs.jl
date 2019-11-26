using JSON
using Formatting
using ArgParse
using Dates
 

function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table s begin

        "--auto-resubmit"
            help = "Whether to give user prompt to resubmitting failed jobs."
            action = :store_true
 
        "--target-year"
            help = "Target year."
            arg_type = Int64
            required = true 

        "--status-file"
            help = "Status file"
            arg_type = String
            default = "JOB_STATUS" 

        "--qsub-options"
            help = "Extra qsub options if needed."
            arg_type = String
            default = ""

    end

    return parse_args(ARGS, s)
end

parsed = parse_commandline()

println(format("# Current real time:  {:s}", Dates.format(now(), "yyyy/mm/dd HH:MM:SS")))

archive_dir = "/glade/scratch/tienyiao/archive"
cases_dir = "./cases"
all_casenames = split(read(pipeline(`ls $cases_dir`), String), "\n", keepempty=false)
working_casenames = split(read(pipeline(`qstat -f`, ignorestatus(`grep Job_Name`), `awk '{print $3}'`), String), "\n", keepempty=false)

println("##### All casenames: #####")
JSON.print(all_casenames, 4)

println("##### Working casenames: #####")
JSON.print(working_casenames, 4)


failed_casenames = []
ok_casenames = []
for (i, casename) in enumerate(all_casenames)

    if casename in working_casenames
        continue
    end

    target_file = joinpath(
        archive_dir,
        casename,
        "ocn", "hist",
        format("{:s}.ocn.h.monthly.{:04d}.nc", casename, parsed["target-year"])
    )

    push!( ( isfile(target_file) ) ? ok_casenames : failed_casenames, casename)
end

if length(failed_casenames) > 0
    println("##### Failed casenames: #####")
    JSON.print(failed_casenames, 4)

    working_dir = pwd()
    println("Do you want to qsub these cases? (yes/no)")
    if parsed["auto-resubmit"] || strip(readline(stdin)) == "yes"
        for casename in failed_casenames
            cd(working_dir)
            println("Submitting ", casename)
            cd(joinpath(cases_dir, casename))
            run(`qsub $(parsed["qsub-options"]) $(casename).run`)
        end
    end
end

open(parsed["status-file"], "w") do io
    write(io, (length(ok_casenames) == length(all_casenames) ? "DONE" : "UNDONE") )
end

if length(ok_casenames) == length(all_casenames)
    println("###### All cases are done  #####")
end
