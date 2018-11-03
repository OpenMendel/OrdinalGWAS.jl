using CSV, PolrGWAS, StatsModels

const datadir = joinpath(dirname(@__FILE__), "..", "data")

@time polrgwas(
@formula(trait ~ 0 + sex),
datadir * "/covariate.txt",
datadir * "/hapmap3",
covtype = [String, String, String, String, Float64, Int],
test = :score
)

@time polrgwas(
@formula(trait ~ 0 + sex),
datadir * "/covariate.txt",
datadir * "/hapmap3",
covtype = [String, String, String, String, Float64, Int],
test = :LRT
)

rm("polrgwas.nullmodel.txt")
rm("polrgwas.scoretest.txt")
rm("polrgwas.lrttest.txt")
