__precompile__()

module OrdinalGWAS

using LinearAlgebra
using CSV, DataFrames, Distributions, Reexport
using SnpArrays, VCFTools, GeneticVariation, BGEN
@reexport using OrdinalMultinomialModels

export ordinalgwas

include("gwas.jl")

end
