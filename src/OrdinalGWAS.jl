__precompile__()

module OrdinalGWAS

using LinearAlgebra
using CSV, DataFrames, Distributions, Reexport, SnpArrays, VCFTools, GeneticVariation
@reexport using OrdinalMultinomialModels

export ordinalgwas

include("gwas.jl")

end
