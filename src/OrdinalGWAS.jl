__precompile__()

module OrdinalGWAS

using LinearAlgebra
using CSV, DataFrames, Distributions, Reexport
using SnpArrays, VCFTools, VariantCallFormat, BGEN
using MathOptInterface
const MOI = MathOptInterface
@reexport using OrdinalMultinomialModels

export ordinalgwas

include("gwas.jl")

end
