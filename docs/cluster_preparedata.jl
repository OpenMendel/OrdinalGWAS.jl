#!/usr/local/bin/julia
#
# This script prepares a data set in current folder. 
# For each of chromosome 1-23, there is a set gzipped Plink files:
# hapmap3.chr.1.bed.gz, hapmap3.chr.1.bim.gz, hapmap3.chr.1.fam.gz
# hapmap3.chr.2.bed.gz, hapmap3.chr.2.bim.gz, hapmap3.chr.2.fam.gz
# ...
# hapmap3.chr.23.bed.gz, hapmap3.chr.23.bim.gz, hapmap3.chr.23.fam.gz
# There is also a csv file "covariate.txt" that contains trait and covariates.
#

# install and load Julia packages
using Pkg
haskey(Pkg.installed(), "SnpArrays") || 
Pkg.add(PackageSpec(url="https://github.com/OpenMendel/SnpArrays.jl.git"))
haskey(Pkg.installed(), "OrdinalMultinomialModels") || 
Pkg.add(PackageSpec(url="https://github.com/OpenMendel/OrdinalMultinomialModels.jl.git"))
haskey(Pkg.installed(), "OrdinalGWAS") || 
Pkg.add(PackageSpec(url="https://github.com/OpenMendel/OrdinalGWAS.jl.git"))
using OrdinalMultinomialModels, OrdinalGWAS, SnpArrays

# split hapmap3 data according to chromosome
datadir = normpath(joinpath(dirname(pathof(OrdinalGWAS)), "../data/"))
SnpArrays.split_plink(datadir * "hapmap3", :chromosome; prefix = "hapmap3.chr.")
# compresse Plink files for chromosome 1-23
for chr in 1:23
    plinkfile = "hapmap3.chr." * string(chr)
    SnpArrays.compress_plink(plinkfile)
end
# delete uncompressed chromosome Plink files
for chr in 1:26
    plinkfile = "hapmap3.chr." * string(chr)
    isfile(plinkfile * ".bed") && rm(plinkfile * ".bed")
    isfile(plinkfile * ".bim") && rm(plinkfile * ".bim")
    isfile(plinkfile * ".fam") && rm(plinkfile * ".fam")
end
# copy covariate.txt file
cp(datadir * "covariate.txt", joinpath(pwd(), "covariate.txt"))
