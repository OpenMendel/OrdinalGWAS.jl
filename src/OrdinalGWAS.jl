__precompile__()

module OrdinalGWAS

using CSV, DataFrames, Distributions, Reexport, SnpArrays
@reexport using OrdinalMultinomialModels

export ordinalgwas

include("gwas.jl")

end
