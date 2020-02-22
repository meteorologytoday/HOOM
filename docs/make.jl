using Documenter

makedocs(
    root    = joinpath(@__DIR__),
    source  = joinpath(@__DIR__, "..", "src"),
    build   = "build",
    clean   = true,
    doctest = true,
    modules = Module[],
    repo    = "",
    highlightsig = true,
    sitename = "HOOM_site",
    expandfirst = [],
)
