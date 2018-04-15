module PolrGWAS

using CSV, DataFrames, PolrModels, SnpArrays, StatsModels

export
    polrgwas

include("polrgwas_score.jl")

end
