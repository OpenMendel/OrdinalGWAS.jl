module PolrGWAS

using CSV, DataFrames, Reexport, SnpArrays
@reexport using PolrModels


export
    polrgwas

include("polrgwas_score.jl")

end
