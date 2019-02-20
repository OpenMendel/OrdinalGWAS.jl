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
if haskey(Pkg.installed(), "SnpArrays")
    Pkg.update("SnpArrays")
else
    Pkg.add(PackageSpec(url="https://github.com/OpenMendel/SnpArrays.jl.git"))
end
if haskey(Pkg.installed(), "OrdinalMultinomialModels")
    Pkg.update("OrdinalMultinomialModels")
else
    Pkg.add(PackageSpec(url="https://github.com/OpenMendel/OrdinalMultinomialModels.jl.git"))
end
if haskey(Pkg.installed(), "OrdinalGWAS")
    Pkg.update("OrdinalGWAS")
else
    Pkg.add(PackageSpec(url="https://github.com/OpenMendel/OrdinalGWAS.jl.git"))
end
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
    rm(plinkfile * ".bed", force=true)
    rm(plinkfile * ".bim", force=true)
    rm(plinkfile * ".fam", force=true)
end
# copy covariate.txt file
cp(datadir * "covariate.txt", joinpath(pwd(), "covariate.txt"))
