using Documenter, PolrGWAS

ENV["DOCUMENTER_DEBUG"] = "true"

makedocs(
    format = :html,
    sitename = "PolrGWAS"
)

deploydocs(
    repo   = "github.com/OpenMendel/PolrGWAS.jl.git",
    target = "build",
    julia  = "1.0",
    deps   = nothing,
    make   = nothing
)
