using Documenter, OrdinalGWAS

ENV["DOCUMENTER_DEBUG"] = "true"

makedocs(
    format = Documenter.HTML(),
    sitename = "OrdinalGWAS",
    modules = [OrdinalGWAS]
)

deploydocs(
    repo   = "github.com/OpenMendel/OrdinalGWAS.jl.git",
    target = "build"
)
