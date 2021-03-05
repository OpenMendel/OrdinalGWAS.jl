# This was run in Julia 1.5.0

using BGEN, OrdinalMultinomialModels, Random, DataFrames, CSV

bgenfile = "bgen_test.bgen"
bgenfile = "/Users/christophergerman/.julia/dev/OrdinalGWAS/data/bgen_test.bgen"
b = Bgen(BGEN.datadir("example.8bits.bgen"))

people = n_samples(b)
nsnps = n_variants(b)
sampleids = samples(b)
variants = parse_variants(b; from_bgen_start=false)
#sample 3 and 15 and 25 to be causal
v3 = minor_allele_dosage!(b, variants[3]; T=Float64, mean_impute=true)
v25 = minor_allele_dosage!(b, variants[25]; T=Float64, mean_impute=true)


Random.seed!(310)
sex = rand(0:1, people)
X = [sex v3 v25]
θ = [-1.5; 0; 1.5]
β = [0.5; 0.5; -0.5]
y = rpolr(X, β, θ, LogitLink())

df = DataFrame(y = y, sex = sex)
CSV.write("bgen_ex.csv", df)

using OrdinalGWAS

b_rsids = rsids(b)
reps1, rem1 = divrem(nsnps, 10) 
repeats = [fill(10, reps1); fill(rem1, 1)]
genenames = ["gene" * string(i) for i in 1:length(repeats .!= 0)]
annotations = vcat(map((gname, reps) -> 
    fill(gname, reps), 
    genenames, repeats)...)
dfannote = DataFrame(gene = annotations, rsids = b_rsids)
CSV.write("bgen_snpsetfile.txt", dfannote, delim = " ", 
    header = false)

