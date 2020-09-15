# This was run in Julia 1.5.0

using GeneticVariation, VCFTools, OrdinalMultinomialModels, Random, DataFrames, CSV

vcffile = "vcf_test.vcf.gz"
D = convert_ds(Float64, vcffile; key="DS", impute=true, center=false, scale=false)

#find two SNPs that are fairly independent to simulate the trait from
using Statistics
toplist = sortperm(vec(var(D, dims=[1])), rev=true)[1:30]
lastcor = 1
for i in toplist, j in toplist 
    curcor = abs(Statistics.cor(D[:, i], D[:, j]))
    if curcor < lastcor
        lastcor = curcor
        println("$i and $j cor is $curcor")
    end
end

#returns 86 and 656 as uncorrelated markers that still have high variance 

people, snps = nsamples(vcffile), nrecords(vcffile)
Random.seed!(310)
sex = rand(0:1, people)
X = [sex D[:, 86] D[:, 656]]
θ = [-1.5; 0; 1.5]
β = [0.5; 0.5; -0.5]
y = rpolr(X, β, θ, LogitLink())

df = DataFrame(y = y, sex = sex)
CSV.write("vcf_example.csv", df)

#ordinalgwas(@formula(y ~ sex), "vcf_example.csv", "test"; geneticformat = "VCF", vcftype=:DS)
#ordinalgwas(@formula(y ~ sex), "vcf_example.csv", "test"; geneticformat = "VCF", vcftype=:DS, snpinds = [86; 656], test=:LRT)