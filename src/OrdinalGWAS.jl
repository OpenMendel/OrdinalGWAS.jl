__precompile__()

module OrdinalGWAS

using LinearAlgebra
using CSV, DataFrames, Distributions, Reexport, SnpArrays
@reexport using OrdinalMultinomialModels

export ordinalgwas, 
ordinalsnpsetgwas,
ordinalgwasGxE

include("gwas.jl")

end
