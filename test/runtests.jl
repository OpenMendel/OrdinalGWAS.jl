using CSV, PolrGWAS, StatsModels

const datadir = joinpath(dirname(@__FILE__), "..", "data")

@time polrgwas(
    @formula(trait ~ 0 + sex),
    datadir * "/covariate.txt",
    datadir * "/hapmap3";
    covartype = [String, String, String, String, Float64, Int]
    )
