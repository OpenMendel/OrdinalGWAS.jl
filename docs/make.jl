using Documenter

ENV["DOCUMENTER_DEBUG"] = "true"

makedocs(
    format = :html,
    sitename = "PolrGWAS"
)

deploydocs(
    repo   = "github.com/OpenMendel/PolrGWAS.jl.git",
    target = "build"
)
