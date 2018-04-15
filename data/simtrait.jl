# code for simulate polr trait from hapmap3 data

using CSV, DataFrames, PolrModels, SnpArrays

# read in Plink fam file
plkfam = CSV.read("hapmap3.fam";
    header=["famid", "perid", "faid", "moid", "sex", "trait"],
    types=[String, String, String, String, Int, Float64],
    delim=' ',
    datarow=1)

# read in Plink bim file
hapmap = SnpArray("hapmap3")
# set SNP 2, 3, 4, 5 as causal snps with effect size βgeno
Xgeno = convert(Matrix{Float64}, hapmap[:, 2:5])
# set sex effect to be 0.5
X = [plkfam[:sex] Xgeno]

# simular ordered categorical trait
srand(123)
J = 4
θ = collect(1.0:J-1)
β = [0.5; 1.25ones(size(Xgeno, 2))]
link = LogitLink()
Y = rpolr(X, β, θ, link)

# write csv file
covdf = [plkfam[1:5] DataFrame(trait=Y)]
CSV.write("covariate.txt", covdf)
