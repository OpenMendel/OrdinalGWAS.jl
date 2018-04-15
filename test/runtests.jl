using CSV, PolrGWAS, StatsModels

const datadir = joinpath(dirname(@__FILE__), "..", "data")

@time polrgwas(
    @formula(trait ~ 0 + sex),
    plinkfile = datadir * "/hapmap3",
    covarfile = datadir * "/covariate.txt",
    covartype = [String, String, String, String, Float64, Int]
    )
