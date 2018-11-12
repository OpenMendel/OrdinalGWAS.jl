using PolrGWAS, StatsModels

const datadir = joinpath(dirname(@__FILE__), "..", "data")

@time polrgwas(
@formula(trait ~ 0 + sex),
datadir * "/covariate.txt",
datadir * "/hapmap3",
test = :score
)

@time polrgwas(
@formula(trait ~ 0 + sex),
datadir * "/covariate.txt",
datadir * "/hapmap3",
test = :LRT
)

rm("polrgwas.nullmodel.txt")
rm("polrgwas.scoretest.txt")
rm("polrgwas.lrttest.txt")
