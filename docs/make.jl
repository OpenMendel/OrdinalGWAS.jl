using Documenter, PolrGWAS

ENV["DOCUMENTER_DEBUG"] = "true"

makedocs(
    format = :html,
    sitename = "PolrGWAS",
    modules = [PolrGWAS]
)

deploydocs(
    repo   = "github.com/OpenMendel/PolrGWAS.jl.git",
    target = "build"
)
