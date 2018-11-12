__precompile__()

module PolrGWAS

using CSV, DataFrames, Distributions, Reexport, SnpArrays
@reexport using PolrModels

export polrgwas

include("gwas.jl")

end
