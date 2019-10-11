# OrdinalGWAS.jl

OrdinalGWAS.jl is a Julia package for performing genome-wide association studies (GWAS) for ordered categorical phenotypes using [proportional odds model](https://en.wikipedia.org/wiki/Ordered_logit) or [ordred Probit model](https://en.wikipedia.org/wiki/Ordered_probit). It is useful when the phenotype takes ordered discrete values, e.g., disease status (undiagnosed, pre-disease, mild, moderate, severe).

## Installation

This package requires Julia v0.7 or later and two other unregistered packages SnpArrays.jl and OrdinalMultinomialModels.jl. The package has not yet been registered and must be installed using the repository location. Start julia and use the ] key to switch to the package manager REPL
```julia
(v1.1) pkg> add https://github.com/OpenMendel/SnpArrays.jl
(v1.1) pkg> add https://github.com/OpenMendel/OrdinalMultinomialModels.jl
(v1.1) pkg> add https://github.com/OpenMendel/OrdinalGWAS.jl
```


```julia
# machine information for this tutorial
versioninfo()
```

    Julia Version 1.2.0
    Commit c6da87ff4b (2019-08-20 00:03 UTC)
    Platform Info:
      OS: macOS (x86_64-apple-darwin18.6.0)
      CPU: Intel(R) Core(TM) i7-4850HQ CPU @ 2.30GHz
      WORD_SIZE: 64
      LIBM: libopenlibm
      LLVM: libLLVM-6.0.1 (ORCJIT, haswell)



```julia
# for use in this tutorial
using BenchmarkTools, CSV, Glob, OrdinalGWAS, SnpArrays
```

## Example data set

`data` folder of the package contains an example data set. In general, user can locate this folder by command


```julia
using OrdinalGWAS
const datadir = normpath(joinpath(dirname(pathof(OrdinalGWAS)), "../data/"))
```




    "/Users/christophergerman/.julia/dev/OrdinalGWAS/data/"




```julia
# content of the data folder
readdir(glob"*.*", datadir)
```




    7-element Array{String,1}:
     "/Users/christophergerman/.julia/dev/OrdinalGWAS/data/covariate.txt"        
     "/Users/christophergerman/.julia/dev/OrdinalGWAS/data/hapmap3.bed"          
     "/Users/christophergerman/.julia/dev/OrdinalGWAS/data/hapmap3.bim"          
     "/Users/christophergerman/.julia/dev/OrdinalGWAS/data/hapmap3.fam"          
     "/Users/christophergerman/.julia/dev/OrdinalGWAS/data/hapmap3.map"          
     "/Users/christophergerman/.julia/dev/OrdinalGWAS/data/hapmap_snpsetfile.txt"
     "/Users/christophergerman/.julia/dev/OrdinalGWAS/data/simtrait.jl"          



## Basic usage

The following command performs GWAS using the [proportional odds model](https://en.wikipedia.org/wiki/Ordered_logit). The output is the fitted null model.


```julia
ordinalgwas(@formula(trait ~ sex), datadir * "covariate.txt", datadir * "hapmap3")
```




    StatsModels.TableRegressionModel{OrdinalMultinomialModel{Int64,Float64,LogitLink},Array{Float64,2}}
    
    trait ~ sex
    
    Coefficients:
    ──────────────────────────────────────────────────────
                   Estimate  Std.Error   t value  Pr(>|t|)
    ──────────────────────────────────────────────────────
    intercept1|2  -1.48564    0.358891  -4.13952    <1e-4 
    intercept2|3  -0.569479   0.341044  -1.66981    0.0959
    intercept3|4   0.429815   0.339642   1.26549    0.2066
    sex            0.424656   0.213914   1.98517    0.0480
    ──────────────────────────────────────────────────────



For documentation of the `ordinalgwas` function, type `?ordinalgwas` in Julia REPL.
```@docs
ordinalgwas
```

### Formula for null model

The first argument specifies the null model without SNP effects, e.g., `@formula(trait ~ sex)`.

### Input files

`ordinalgwas` expects two input files: one for responses plus covariates (second argument), the other the Plink files for genotypes (third argument).

#### Covariate and trait file

Covariates and phenotype are provided in a csv file, e.g., `covariate.txt`, which has one header line for variable names. In this example, variable `trait` is the ordered categorical phenotypes coded as integers 1 to 4. We want to include variable `sex` as the covariate in GWAS.


```julia
run(`head $(datadir)covariate.txt`);
```

    famid,perid,faid,moid,sex,trait
    2431,NA19916,0,0,1,4
    2424,NA19835,0,0,2,4
    2469,NA20282,0,0,2,4
    2368,NA19703,0,0,1,3
    2425,NA19901,0,0,2,3
    2427,NA19908,0,0,1,4
    2430,NA19914,0,0,2,4
    2470,NA20287,0,0,2,1
    2436,NA19713,0,0,2,3


#### Plink file

Genotype data is available as binary Plink files.


```julia
readdir(glob"hapmap3.*", datadir)
```




    4-element Array{String,1}:
     "/Users/christophergerman/.julia/dev/OrdinalGWAS/data/hapmap3.bed"
     "/Users/christophergerman/.julia/dev/OrdinalGWAS/data/hapmap3.bim"
     "/Users/christophergerman/.julia/dev/OrdinalGWAS/data/hapmap3.fam"
     "/Users/christophergerman/.julia/dev/OrdinalGWAS/data/hapmap3.map"



In this example, there are 324 samples at 13,928 SNPs.


```julia
size(SnpArray(datadir * "hapmap3.bed"))
```




    (324, 13928)



Compressed Plink files are supported. For example, if Plink files are `hapmap3.bed.gz`, `hapmap3.bim.gz` and `hapmap3.fam.gz`, the same command
```julia
ordinalgwas(@formula(trait ~ sex), datadir * "covariate.txt", datadir * "hapmap3")
```
still works. Check all supported compression format by


```julia
SnpArrays.ALLOWED_FORMAT
```




    6-element Array{String,1}:
     "gz"  
     "zlib"
     "zz"  
     "xz"  
     "zst" 
     "bz2" 



### Output files

`ordinalgwas` outputs two files: `ordinalgwas.null.txt` and `ordinalgwas.pval.txt`. 

* `ordinalgwas.null.txt` lists the estimated null model (without SNPs). 


```julia
run(`cat ordinalgwas.null.txt`);
```

    StatsModels.TableRegressionModel{OrdinalMultinomialModel{Int64,Float64,LogitLink},Array{Float64,2}}
    
    trait ~ sex
    
    Coefficients:
    ──────────────────────────────────────────────────────
                   Estimate  Std.Error   t value  Pr(>|t|)
    ──────────────────────────────────────────────────────
    intercept1|2  -1.48564    0.358891  -4.13952    <1e-4 
    intercept2|3  -0.569479   0.341044  -1.66981    0.0959
    intercept3|4   0.429815   0.339642   1.26549    0.2066
    sex            0.424656   0.213914   1.98517    0.0480
    ──────────────────────────────────────────────────────

* `ordinalgwas.pval.txt` tallies the SNPs and their pvalues. 


```julia
run(`head ordinalgwas.pval.txt`);
```

    chr,pos,snpid,maf,hwepval,pval
    1,554484,rs10458597,0.0,1.0,1.0
    1,758311,rs12562034,0.07763975155279501,0.40987633326666817,0.004565312839540994
    1,967643,rs2710875,0.32407407407407407,4.0762491007057454e-7,3.108283828554874e-5
    1,1168108,rs11260566,0.19158878504672894,0.12856822794469086,1.2168672367668889e-5
    1,1375074,rs1312568,0.441358024691358,2.537601965061486e-19,0.008206860046175225
    1,1588771,rs35154105,0.0,1.0,1.0
    1,1789051,rs16824508,0.00462962962962965,0.9332783156468178,0.5111981332544
    1,1990452,rs2678939,0.4537037037037037,5.0769595770843126e-11,0.2997282957184783
    1,2194615,rs7553178,0.22685185185185186,0.17056143157457782,0.1713331245805063


Output file names can be changed by the `nullfile` and `pvalfile` keywords respectively. For example, 
```julia
ordinalgwas(@formula(trait ~ sex), datadir * "covariate.txt", datadir * "hapmap3", pvalfile="ordinalgwas.pval.txt.gz")
```
will output the p-value file in compressed gz format.

### Subsamples

Use the keyword `covrowinds` to specify selected samples in the covarite file. Use the keyword `bedrowinds` to specify selected samples in the Plink bed file. For example, to use the first 300 samples in both covariate and bed file:
```julia
ordinalgwas(@formula(trait ~ sex), covfile, plkfile, covrowinds=1:300, bedrowinds=1:300)
```
!!! note

    Users should always make sure that the selected samples in covariate file match exactly those in bed file. 

### Input non-genetic data as DataFrame

Internally `ordinalgwas` parses the covariate file as a DataFrame by `CSV.read(covfile)`. For covariate file of other formats, users can parse it as a DataFrame and then input the DataFrame to `ordinalgwas` directly.
```julia
ordinalgwas(@formula(trait ~ sex), df, plinkfile)
```
!!! note

    Users should always make sure that individuals in covariate file or DataFrame match those in Plink fam file. 

For example, following code checks that the first 2 columns of the `covariate.txt` file match the first 2 columns of the `hapmap3.fam` file exactly.


```julia
covdf = CSV.read(datadir * "covariate.txt")
plkfam = CSV.read(datadir * "hapmap3.fam", header=0, delim=' ')
all(covdf[!, 1] .== plkfam[!, 1]) && all(covdf[!, 2] .== plkfam[!, 2])
```




    true



### Timing

For this moderate-sized data set, `ordinalgwas` takes less than 0.2 second.


```julia
@btime(ordinalgwas(@formula(trait ~ sex), datadir * "covariate.txt", datadir * "hapmap3"));
```

      201.963 ms (639651 allocations: 32.71 MiB)



```julia
# clean up
rm("ordinalgwas.null.txt", force=true)
rm("ordinalgwas.pval.txt", force=true)
```

## Link functions

The `link` keyword argument of `ordinalgwas` can take value:  
- `LogitLink()`, proportional odds model (default),  
- `ProbitLink()`, ordred Probit model,  
- `CloglogLink()`, proportional hazards model, or 
- `CauchyLink()`.

For example, to perform GWAS using the ordred Probit model


```julia
ordinalgwas(@formula(trait ~ sex), datadir * "covariate.txt", datadir * "hapmap3", 
    link=ProbitLink(), nullfile="opm.null.txt", pvalfile="opm.pval.txt")
```




    StatsModels.TableRegressionModel{OrdinalMultinomialModel{Int64,Float64,ProbitLink},Array{Float64,2}}
    
    trait ~ sex
    
    Coefficients:
    ──────────────────────────────────────────────────────
                   Estimate  Std.Error   t value  Pr(>|t|)
    ──────────────────────────────────────────────────────
    intercept1|2  -0.866156   0.210677  -4.11129    <1e-4 
    intercept2|3  -0.359878   0.205817  -1.74854    0.0813
    intercept3|4   0.247054   0.205382   1.2029     0.2299
    sex            0.251058   0.128225   1.95795    0.0511
    ──────────────────────────────────────────────────────



The estimates in null model and p-values are slightly different from those in proportional odds moodel.


```julia
run(`head opm.pval.txt`);
```

    chr,pos,snpid,maf,hwepval,pval
    1,554484,rs10458597,0.0,1.0,1.0
    1,758311,rs12562034,0.07763975155279501,0.40987633326666817,0.010076916742300138
    1,967643,rs2710875,0.32407407407407407,4.0762491007057454e-7,2.6272564941853933e-5
    1,1168108,rs11260566,0.19158878504672894,0.12856822794469086,1.089748485107844e-5
    1,1375074,rs1312568,0.441358024691358,2.537601965061486e-19,0.005102883990438149
    1,1588771,rs35154105,0.0,1.0,1.0
    1,1789051,rs16824508,0.00462962962962965,0.9332783156468178,0.48653776297859236
    1,1990452,rs2678939,0.4537037037037037,5.0769595770843126e-11,0.33231290090455434
    1,2194615,rs7553178,0.22685185185185186,0.17056143157457782,0.2591551397719743



```julia
rm("opm.null.txt", force=true)
rm("opm.pval.txt", force=true)
```

## SNP models

Genotypes are translated into numeric values according to different genetic model, which is specified by the `snpmodel` keyword. Default is `ADDITIVE_MODEL`.

| Genotype | `SnpArray` | `ADDITIVE_MODEL` | `DOMINANT_MODEL` | `RECESSIVE_MODEL` |    
|:---:|:---:|:---:|:---:|:---:|  
| A1,A1 | 0x00 | 0 | 0 | 0 |  
| missing | 0x01 | NaN | NaN | NaN |
| A1,A2 | 0x02 | 1 | 1 | 0 |  
| A2,A2 | 0x03 | 2 | 1 | 1 |  

!!! note

    `ordinalgwas` imputes missing genotypes according to minor allele frequencies. 
    
Users are advised to impute genotypes using more sophiscated methods before GWAS.

## SNP and/or sample masks

In practice, we often perform GWAS on selected SNPs and/or selected samples. They can be specified by the `snpinds`, `covrowinds` and `bedrowinds` keywords of `ordinalgwas` function. 

For example, to perform GWAS on SNPs with minor allele frequency (MAF) above 0.05


```julia
# create SNP mask
snpinds = maf(SnpArray("../data/hapmap3.bed")) .≥ 0.05
# GWAS on selected SNPs
@time ordinalgwas(@formula(trait ~ sex), datadir * "covariate.txt", datadir * "hapmap3", 
    snpinds=snpinds, nullfile="commonvariant.null.txt", pvalfile="commonvariant.pval.txt")
```

      0.202024 seconds (562.98 k allocations: 28.824 MiB, 5.82% gc time)





    StatsModels.TableRegressionModel{OrdinalMultinomialModel{Int64,Float64,LogitLink},Array{Float64,2}}
    
    trait ~ sex
    
    Coefficients:
    ──────────────────────────────────────────────────────
                   Estimate  Std.Error   t value  Pr(>|t|)
    ──────────────────────────────────────────────────────
    intercept1|2  -1.48564    0.358891  -4.13952    <1e-4 
    intercept2|3  -0.569479   0.341044  -1.66981    0.0959
    intercept3|4   0.429815   0.339642   1.26549    0.2066
    sex            0.424656   0.213914   1.98517    0.0480
    ──────────────────────────────────────────────────────




```julia
run(`head commonvariant.pval.txt`);
```

    chr,pos,snpid,maf,hwepval,pval
    1,758311,rs12562034,0.07763975155279501,0.40987633326666817,0.004565312839540994
    1,967643,rs2710875,0.32407407407407407,4.0762491007057454e-7,3.108283828554874e-5
    1,1168108,rs11260566,0.19158878504672894,0.12856822794469086,1.2168672367668889e-5
    1,1375074,rs1312568,0.441358024691358,2.537601965061486e-19,0.008206860046175225
    1,1990452,rs2678939,0.4537037037037037,5.0769595770843126e-11,0.2997282957184783
    1,2194615,rs7553178,0.22685185185185186,0.17056143157457782,0.1713331245805063
    1,2396747,rs13376356,0.1448598130841121,0.905307921507814,0.5320416198875456
    1,2823603,rs1563468,0.4830246913580247,4.2306553724392545e-9,0.225191391783573
    1,3025087,rs6690373,0.2538699690402477,9.23864188719278e-8,0.7018469417717486



```julia
# extra header line in commonvariant.pval.txt
countlines("commonvariant.pval.txt"), count(snpinds)
```




    (12086, 12085)




```julia
# clean up
rm("commonvariant.null.txt", force=true)
rm("commonvariant.pval.txt", force=true)
```

`covrowinds` specify the samples in the covariate file and `bedrowinds` for SnpArray. User should be particularly careful when these two keyword. Selected rows in SnpArray should exactly match the samples in the null model. Otherwise the results are meaningless.

## Likelihood ratio test (LRT)

By default, `ordinalgwas` calculates p-value for each SNP using score test. Score test is fast because it doesn't require fitting alternative model for each SNP. User can request likelihood ratio test (LRT) using keyword `test=:lrt`. LRT is much slower but may be more powerful than score test.


```julia
@time ordinalgwas(@formula(trait ~ sex), datadir * "covariate.txt", datadir * "hapmap3", 
    test=:LRT, nullfile="lrt.null.txt", pvalfile="lrt.pval.txt")
```

     27.804879 seconds (8.59 M allocations: 2.106 GiB, 1.77% gc time)





    StatsModels.TableRegressionModel{OrdinalMultinomialModel{Int64,Float64,LogitLink},Array{Float64,2}}
    
    trait ~ sex
    
    Coefficients:
    ──────────────────────────────────────────────────────
                   Estimate  Std.Error   t value  Pr(>|t|)
    ──────────────────────────────────────────────────────
    intercept1|2  -1.48564    0.358891  -4.13952    <1e-4 
    intercept2|3  -0.569479   0.341044  -1.66981    0.0959
    intercept3|4   0.429815   0.339642   1.26549    0.2066
    sex            0.424656   0.213914   1.98517    0.0480
    ──────────────────────────────────────────────────────



Note the extra `effect` column in pvalfile, which is the effect size (regression coefficient) for each SNP. 


```julia
run(`head lrt.pval.txt`);
```

    chr,pos,snpid,maf,hwepval,effect,pval
    1,554484,rs10458597,0.0,1.0,0.0,1.0
    1,758311,rs12562034,0.07763975155279501,0.40987633326666817,-1.0057833719544331,0.0019185836579804134
    1,967643,rs2710875,0.32407407407407407,4.0762491007057454e-7,-0.6488560566295055,1.805050556976241e-5
    1,1168108,rs11260566,0.19158878504672894,0.12856822794469086,-0.9157225669357879,5.873384712685568e-6
    1,1375074,rs1312568,0.441358024691358,2.537601965061486e-19,-0.33181366525772593,0.008081022577832324
    1,1588771,rs35154105,0.0,1.0,0.0,1.0
    1,1789051,rs16824508,0.00462962962962965,0.9332783156468178,-0.7338026388701573,0.5169027130129711
    1,1990452,rs2678939,0.4537037037037037,5.0769595770843126e-11,-0.13586499231819726,0.29946402200912603
    1,2194615,rs7553178,0.22685185185185186,0.17056143157457782,-0.2512075640440123,0.16151069094439868



```julia
# clean up
rm("lrt.pval.txt", force=true)
rm("lrt.null.txt", force=true)
```

In this example, GWAS by score test takes less than 0.2 second, while GWAS by LRT takes about 20 seconds. About 100 fold difference in run time. 

## Score test for screening, LRT for power 

For large data sets, a practical solution is to perform score test first, then re-do LRT for the most promising SNPs according to score test p-values.

**Step 1**: Perform score test GWAS, results in `score.pval.txt`.


```julia
@time ordinalgwas(@formula(trait ~ sex), datadir * "covariate.txt", datadir * "hapmap3", 
    test=:score, pvalfile="score.pval.txt");
```

      0.272130 seconds (668.77 k allocations: 34.251 MiB, 10.04% gc time)



```julia
run(`head score.pval.txt`);
```

    chr,pos,snpid,maf,hwepval,pval
    1,554484,rs10458597,0.0,1.0,1.0
    1,758311,rs12562034,0.07763975155279501,0.40987633326666817,0.004565312839540994
    1,967643,rs2710875,0.32407407407407407,4.0762491007057454e-7,3.108283828554874e-5
    1,1168108,rs11260566,0.19158878504672894,0.12856822794469086,1.2168672367668889e-5
    1,1375074,rs1312568,0.441358024691358,2.537601965061486e-19,0.008206860046175225
    1,1588771,rs35154105,0.0,1.0,1.0
    1,1789051,rs16824508,0.00462962962962965,0.9332783156468178,0.5111981332544
    1,1990452,rs2678939,0.4537037037037037,5.0769595770843126e-11,0.2997282957184783
    1,2194615,rs7553178,0.22685185185185186,0.17056143157457782,0.1713331245805063


**Step 2**: Sort score test p-values and find top 10 SNPs.


```julia
scorepvals = CSV.read("score.pval.txt")[!, 5] # p-values in 5th column
tophits = sortperm(scorepvals)[1:10] # indices of 10 SNPs with smallest p-values
scorepvals[tophits] # smallest 10 p-values
```




    10-element Array{Float64,1}:
     1.9481897837873706e-72
     1.9481897837874264e-72
     1.9481897837874264e-72
     1.9481897837874264e-72
     1.948189783787482e-72 
     1.948189783787482e-72 
     3.216959988402954e-72 
     3.216959988402954e-72 
     3.216959988403046e-72 
     3.216959988403046e-72 



**Step 3**: Re-do LRT on top hits.


```julia
@time ordinalgwas(@formula(trait ~ sex), datadir * "covariate.txt", datadir * "hapmap3", 
    snpinds=tophits, test=:LRT, pvalfile="lrt.pval.txt");
```

      0.036746 seconds (61.81 k allocations: 3.704 MiB)



```julia
run(`cat lrt.pval.txt`);
```

    chr,pos,snpid,maf,hwepval,effect,pval
    2,111915648,rs1004117,0.009287925696594423,3.216959988402954e-72,0.29754964794084027,0.5418602577787199
    9,45757243,rs2772909,0.0092592592592593,1.9481897837873706e-72,1.2617886438731196,0.06079130674507882
    9,102275751,rs17853384,0.003095975232198178,3.216959988403046e-72,0.4668871230700536,0.5300792023608788
    16,28673672,rs3874769,0.003095975232198178,3.216959988403046e-72,0.46650749274735515,0.530436201884302
    16,44943958,rs4640368,0.006172839506172867,1.9481897837874264e-72,0.17634494923727342,0.7866349302676543
    17,23519525,rs34876265,0.0030864197530864335,1.948189783787482e-72,-5.2012149325995,0.3138904865441721
    19,42449922,rs2918343,0.009287925696594423,3.216959988402954e-72,0.6394911477064131,0.15203146280636515
    21,43185761,rs8130793,0.006172839506172867,1.9481897837874264e-72,-4.636834541616616,0.1536110542875311
    23,137700451,rs7061097,0.0030864197530864335,1.948189783787482e-72,0.24982571474992607,0.7343368996308214
    26,410,rs28412942,0.006172839506172867,1.9481897837874264e-72,-4.6368345416179055,0.1536110542875311



```julia
# clean up
rm("ordinalgwas.null.txt", force=true)
rm("score.pval.txt", force=true)
rm("lrt.pval.txt", force=true)
```

## GxE or other interactions

In many applications, we want to test SNP effect and/or its interaction with other terms. `testformula` keyword specifies the test unit **besides** the covariates in `nullformula`. 

In following example, keyword `testformula=@formula(trait ~ snp + snp & sex)` instructs `ordinalgwas` to test joint effect of `snp` and `snp & sex` interaction.


```julia
ordinalgwas(@formula(trait ~ sex), datadir * "covariate.txt", datadir * "hapmap3", 
    pvalfile="GxE.pval.txt", testformula=@formula(trait ~ snp + snp & sex));
```


```julia
run(`head GxE.pval.txt`);
```

    chr,pos,snpid,maf,hwepval,pval
    1,554484,rs10458597,0.0,1.0,1.0
    1,758311,rs12562034,0.07763975155279501,0.40987633326666817,0.01744601041225423
    1,967643,rs2710875,0.32407407407407407,4.0762491007057454e-7,0.0001667073239489097
    1,1168108,rs11260566,0.19158878504672894,0.12856822794469086,4.763762457893383e-5
    1,1375074,rs1312568,0.441358024691358,2.537601965061486e-19,0.02913847124299368
    1,1588771,rs35154105,0.0,1.0,1.0
    1,1789051,rs16824508,0.00462962962962965,0.9332783156468178,0.29643631149443217
    1,1990452,rs2678939,0.4537037037037037,5.0769595770843126e-11,0.37924580479348796
    1,2194615,rs7553178,0.22685185185185186,0.17056143157457782,0.325582269932396



```julia
# clean up
rm("ordinalgwas.null.txt")
rm("GxE.pval.txt")
```

For some applications, the user may want to simply test the GxE interaction effeect. This requires fitting the snp in the null model and is quite slower, but the command `ordinalgwasGxE()` can be used test the interaction effect.
To do this you must specify the environmental variable. 

For documentation of the `ordinalgwasGxE` function, type `?ordinalgwasGxE` in Julia REPL.
```@docs
ordinalgwasGxE
```


```julia
ordinalgwasGxE(@formula(trait ~ sex), datadir * "covariate.txt", datadir * "hapmap3",
    :sex, pvalfile = "gxe_snp.pval.txt", snpinds=1:5, test=:score)
```




    StatsModels.TableRegressionModel{OrdinalMultinomialModel{Int64,Float64,LogitLink},Array{Float64,2}}
    
    trait ~ sex
    
    Coefficients:
    ──────────────────────────────────────────────────────
                   Estimate  Std.Error   t value  Pr(>|t|)
    ──────────────────────────────────────────────────────
    intercept1|2  -1.48564    0.358891  -4.13952    <1e-4 
    intercept2|3  -0.569479   0.341044  -1.66981    0.0959
    intercept3|4   0.429815   0.339642   1.26549    0.2066
    sex            0.424656   0.213914   1.98517    0.0480
    ──────────────────────────────────────────────────────




```julia
run(`head gxe_snp.pval.txt`);
```

    chr,pos,snpid,maf,hwepval,snpeffectnull,pval
    1,554484,rs10458597,0.0,1.0,0.0,1.0
    1,758311,rs12562034,0.07763975155279501,0.40987633326666817,-1.0057833719544331,0.6377422425977491
    1,967643,rs2710875,0.32407407407407407,4.0762491007057454e-7,-0.6488560566295055,0.9667114198051628
    1,1168108,rs11260566,0.19158878504672894,0.12856822794469086,-0.9157225669357879,0.26352674694121003
    1,1375074,rs1312568,0.441358024691358,2.537601965061486e-19,-0.33181366525772593,0.7811133315582837



```julia
# clean up
rm("gxe_snp.pval.txt", force=true)
```

## SNP-set testing

In many applications, we want to test SNP-set. The function `ordinalsnpsetgwas()` can be used to do this. The snpset can be specified as either a window (test every X snps), a filename that specifies an input file, with no header and two columns separated by a space. The first column must contain the snpset ID and the second column must contain the snpid's identical to the bimfile, or an AbstractVector that allows you to specify the snps you want to perform one joint snpset test for.

For documentation of the `ordinalsnpsetgwas` function, type `?ordinalsnpsetgwas` in Julia REPL.
```@docs
ordinalsnpsetgwas
```

In following example, we perform a snpset test on the 50th to 55th snps. 


```julia
ordinalsnpsetgwas(@formula(trait ~ sex), datadir * "covariate.txt", datadir * "hapmap3",
    pvalfile = "snpset.pval.txt", snpset = 50:55)
```




    StatsModels.TableRegressionModel{OrdinalMultinomialModel{Int64,Float64,LogitLink},Array{Float64,2}}
    
    trait ~ sex
    
    Coefficients:
    ──────────────────────────────────────────────────────
                   Estimate  Std.Error   t value  Pr(>|t|)
    ──────────────────────────────────────────────────────
    intercept1|2  -1.48564    0.358891  -4.13952    <1e-4 
    intercept2|3  -0.569479   0.341044  -1.66981    0.0959
    intercept3|4   0.429815   0.339642   1.26549    0.2066
    sex            0.424656   0.213914   1.98517    0.0480
    ──────────────────────────────────────────────────────




```julia
run(`head snpset.pval.txt`);
```

    The joint pvalue of snps indexed at 50:55 is 0.36471265366639605



```julia
# clean up
rm("snpset.pval.txt", force=true)
rm("ordinalgwas.null.txt", force=true)
```

In the following example we run a snpset test on the annotated snpset file.


```julia
ordinalsnpsetgwas(@formula(trait ~ sex), datadir * "covariate.txt", datadir * "hapmap3",
    pvalfile = "snpset.pval.txt", snpset = datadir * "/hapmap_snpsetfile.txt")
```




    StatsModels.TableRegressionModel{OrdinalMultinomialModel{Int64,Float64,LogitLink},Array{Float64,2}}
    
    trait ~ sex
    
    Coefficients:
    ──────────────────────────────────────────────────────
                   Estimate  Std.Error   t value  Pr(>|t|)
    ──────────────────────────────────────────────────────
    intercept1|2  -1.48564    0.358891  -4.13952    <1e-4 
    intercept2|3  -0.569479   0.341044  -1.66981    0.0959
    intercept3|4   0.429815   0.339642   1.26549    0.2066
    sex            0.424656   0.213914   1.98517    0.0480
    ──────────────────────────────────────────────────────




```julia
run(`head snpset.pval.txt`);
```

    snpsetid,nsnps,pval
    gene1,93,1.72133946991466e-5
    gene2,93,0.03692497684330217
    gene3,93,0.7478549371392197
    gene4,92,0.027650798223816195
    gene5,93,0.6119581594582347
    gene6,93,0.029184642230868463
    gene7,93,0.3916007348891039
    gene8,93,0.1253957410986018
    gene9,93,0.7085635621609013



```julia
# clean up
rm("snpset.pval.txt", force=true)
rm("ordinalgwas.null.txt", force=true)
```

In the following example we run a snpset test on every 15 snps.


```julia
ordinalsnpsetgwas(@formula(trait ~ sex), datadir * "covariate.txt", datadir * "hapmap3",
    pvalfile = "snpset.pval.txt", snpset=15)
```




    StatsModels.TableRegressionModel{OrdinalMultinomialModel{Int64,Float64,LogitLink},Array{Float64,2}}
    
    trait ~ sex
    
    Coefficients:
    ──────────────────────────────────────────────────────
                   Estimate  Std.Error   t value  Pr(>|t|)
    ──────────────────────────────────────────────────────
    intercept1|2  -1.48564    0.358891  -4.13952    <1e-4 
    intercept2|3  -0.569479   0.341044  -1.66981    0.0959
    intercept3|4   0.429815   0.339642   1.26549    0.2066
    sex            0.424656   0.213914   1.98517    0.0480
    ──────────────────────────────────────────────────────




```julia
run(`head snpset.pval.txt`);
```

    startchr,startpos,startsnpid,endchr,endpos,endsnpid,pval
    1,554484,rs10458597,1,3431124,rs12093117,1.9394189465782687e-13
    1,3633945,rs10910017,1,6514524,rs932112,0.11077869162820388
    1,6715827,rs441515,1,9534606,rs4926480,0.2742450817971269
    1,9737551,rs12047054,1,12559747,rs4845907,0.49346113650468615
    1,12760427,rs848577,1,16021797,rs6679870,0.15447358658305502
    1,16228774,rs1972359,1,19100349,rs9426794,0.1544232917025094
    1,19301516,rs4912046,1,22122176,rs9426785,0.474987334946532
    1,22323074,rs2235529,1,25166528,rs4648890,0.41246096621475836
    1,25368553,rs7527379,1,28208168,rs12140070,0.1645813527877886



```julia
# clean up
rm("snpset.pval.txt", force=true)
rm("ordinalgwas.null.txt", force=true)
```

## Plotting Results

To plot the GWAS results, use the [MendelPlots.jl package](https://openmendel.github.io/MendelPlots.jl/latest/).

## Docker

For ease of using OrdinalGWAS, we provide a Dockerfile so users don't need to install Julia and required packages. Only Docker app needs to be installed in order to run analysis. Following is tested on Docker 2.0.0.0-mac78.

**Step 1**: Create a Dockerfile with content [here](https://raw.githubusercontent.com/OpenMendel/OrdinalGWAS.jl/master/docker/Dockerfile), or, if the bash command `wget` is available, obtain Dockerfile by
```bash
# on command line
wget https://raw.githubusercontent.com/OpenMendel/OrdinalGWAS.jl/master/docker/Dockerfile
```

**Step 2**: Build a docker image called `ordinalgwas-app`, assuming that the Dockerfile is located in the `../docker` folder. Building the image for the first time can take up to 10 minutes; but it only needs to be done once.
```bash
# on command line
docker build -t ordinalgwas-app ../docker/
```

**Step 3**: Suppose data files are located at `/path/to/data` folder, run analysis by
```bash
# on command line
docker run -v /path/to/data:/data -t ordinalgwas-app julia -e 'using OrdinalGWAS; ordinalgwas(@formula(trait ~ sex), "/data/covariate.txt", "/data/hapmap3", nullfile="/data/ordinalgwas.null.txt", pvalfile="/data/ordinalgwas");'
```

Here  
- `-t ordinalgwas-app` creates a container using the `ordinalgwas-app` image build in step 2.  
- `-v /path/to/data:/data` maps the `/path/to/data` folder on host machine to the `/data` folder within the container. 
- `julia -e 'using OrdinalGWAS; ordinalgwas(@formula(trait ~ sex), "/data/covariate.txt", "/data/hapmap3", nullfile="/data/ordinalgwas.null.txt", pvalfile="/data/ordinalgwas");` calls Julia and runs `ordinalgwas` function. 

The output files are written in `/path/to/data` directory.

## Multiple Plink file sets

In large scale studies, genotypes data are split into multiple Plink files, e.g., by chromosome. Then GWAS analysis can be done in parallel. This can be achieved by two steps.

Let's first create demo data by splitting hapmap3 according to chromosome:


```julia
# split example hapmap3 data according to chromosome
SnpArrays.split_plink(datadir * "hapmap3", :chromosome; prefix=datadir * "hapmap3.chr.")
readdir(glob"hapmap3.chr.*", datadir)
```




    75-element Array{String,1}:
     "/Users/christophergerman/.julia/dev/OrdinalGWAS/data/hapmap3.chr.1.bed" 
     "/Users/christophergerman/.julia/dev/OrdinalGWAS/data/hapmap3.chr.1.bim" 
     "/Users/christophergerman/.julia/dev/OrdinalGWAS/data/hapmap3.chr.1.fam" 
     "/Users/christophergerman/.julia/dev/OrdinalGWAS/data/hapmap3.chr.10.bed"
     "/Users/christophergerman/.julia/dev/OrdinalGWAS/data/hapmap3.chr.10.bim"
     "/Users/christophergerman/.julia/dev/OrdinalGWAS/data/hapmap3.chr.10.fam"
     "/Users/christophergerman/.julia/dev/OrdinalGWAS/data/hapmap3.chr.11.bed"
     "/Users/christophergerman/.julia/dev/OrdinalGWAS/data/hapmap3.chr.11.bim"
     "/Users/christophergerman/.julia/dev/OrdinalGWAS/data/hapmap3.chr.11.fam"
     "/Users/christophergerman/.julia/dev/OrdinalGWAS/data/hapmap3.chr.12.bed"
     "/Users/christophergerman/.julia/dev/OrdinalGWAS/data/hapmap3.chr.12.bim"
     "/Users/christophergerman/.julia/dev/OrdinalGWAS/data/hapmap3.chr.12.fam"
     "/Users/christophergerman/.julia/dev/OrdinalGWAS/data/hapmap3.chr.13.bed"
     ⋮                                                                        
     "/Users/christophergerman/.julia/dev/OrdinalGWAS/data/hapmap3.chr.6.bed" 
     "/Users/christophergerman/.julia/dev/OrdinalGWAS/data/hapmap3.chr.6.bim" 
     "/Users/christophergerman/.julia/dev/OrdinalGWAS/data/hapmap3.chr.6.fam" 
     "/Users/christophergerman/.julia/dev/OrdinalGWAS/data/hapmap3.chr.7.bed" 
     "/Users/christophergerman/.julia/dev/OrdinalGWAS/data/hapmap3.chr.7.bim" 
     "/Users/christophergerman/.julia/dev/OrdinalGWAS/data/hapmap3.chr.7.fam" 
     "/Users/christophergerman/.julia/dev/OrdinalGWAS/data/hapmap3.chr.8.bed" 
     "/Users/christophergerman/.julia/dev/OrdinalGWAS/data/hapmap3.chr.8.bim" 
     "/Users/christophergerman/.julia/dev/OrdinalGWAS/data/hapmap3.chr.8.fam" 
     "/Users/christophergerman/.julia/dev/OrdinalGWAS/data/hapmap3.chr.9.bed" 
     "/Users/christophergerman/.julia/dev/OrdinalGWAS/data/hapmap3.chr.9.bim" 
     "/Users/christophergerman/.julia/dev/OrdinalGWAS/data/hapmap3.chr.9.fam" 



Step 1: Fit the null model. Setting third argument `plinkfile` to `nothing` instructs `ordinalgwas` function to fit the null model only.


```julia
nm = ordinalgwas(@formula(trait ~ sex), datadir * "covariate.txt", nothing)
```




    StatsModels.TableRegressionModel{OrdinalMultinomialModel{Int64,Float64,LogitLink},Array{Float64,2}}
    
    trait ~ sex
    
    Coefficients:
    ──────────────────────────────────────────────────────
                   Estimate  Std.Error   t value  Pr(>|t|)
    ──────────────────────────────────────────────────────
    intercept1|2  -1.48564    0.358891  -4.13952    <1e-4 
    intercept2|3  -0.569479   0.341044  -1.66981    0.0959
    intercept3|4   0.429815   0.339642   1.26549    0.2066
    sex            0.424656   0.213914   1.98517    0.0480
    ──────────────────────────────────────────────────────



Step 2: GWAS for each chromosome.


```julia
# this part can be submitted as separate jobs
for chr in 1:23
    plinkfile = datadir * "hapmap3.chr." * string(chr)
    pvalfile = plinkfile * ".pval.txt" 
    ordinalgwas(nm, plinkfile, pvalfile=pvalfile)
end
```


```julia
# show the result files
readdir(glob"*.pval.txt", datadir)
```




    23-element Array{String,1}:
     "/Users/christophergerman/.julia/dev/OrdinalGWAS/data/hapmap3.chr.1.pval.txt" 
     "/Users/christophergerman/.julia/dev/OrdinalGWAS/data/hapmap3.chr.10.pval.txt"
     "/Users/christophergerman/.julia/dev/OrdinalGWAS/data/hapmap3.chr.11.pval.txt"
     "/Users/christophergerman/.julia/dev/OrdinalGWAS/data/hapmap3.chr.12.pval.txt"
     "/Users/christophergerman/.julia/dev/OrdinalGWAS/data/hapmap3.chr.13.pval.txt"
     "/Users/christophergerman/.julia/dev/OrdinalGWAS/data/hapmap3.chr.14.pval.txt"
     "/Users/christophergerman/.julia/dev/OrdinalGWAS/data/hapmap3.chr.15.pval.txt"
     "/Users/christophergerman/.julia/dev/OrdinalGWAS/data/hapmap3.chr.16.pval.txt"
     "/Users/christophergerman/.julia/dev/OrdinalGWAS/data/hapmap3.chr.17.pval.txt"
     "/Users/christophergerman/.julia/dev/OrdinalGWAS/data/hapmap3.chr.18.pval.txt"
     "/Users/christophergerman/.julia/dev/OrdinalGWAS/data/hapmap3.chr.19.pval.txt"
     "/Users/christophergerman/.julia/dev/OrdinalGWAS/data/hapmap3.chr.2.pval.txt" 
     "/Users/christophergerman/.julia/dev/OrdinalGWAS/data/hapmap3.chr.20.pval.txt"
     "/Users/christophergerman/.julia/dev/OrdinalGWAS/data/hapmap3.chr.21.pval.txt"
     "/Users/christophergerman/.julia/dev/OrdinalGWAS/data/hapmap3.chr.22.pval.txt"
     "/Users/christophergerman/.julia/dev/OrdinalGWAS/data/hapmap3.chr.23.pval.txt"
     "/Users/christophergerman/.julia/dev/OrdinalGWAS/data/hapmap3.chr.3.pval.txt" 
     "/Users/christophergerman/.julia/dev/OrdinalGWAS/data/hapmap3.chr.4.pval.txt" 
     "/Users/christophergerman/.julia/dev/OrdinalGWAS/data/hapmap3.chr.5.pval.txt" 
     "/Users/christophergerman/.julia/dev/OrdinalGWAS/data/hapmap3.chr.6.pval.txt" 
     "/Users/christophergerman/.julia/dev/OrdinalGWAS/data/hapmap3.chr.7.pval.txt" 
     "/Users/christophergerman/.julia/dev/OrdinalGWAS/data/hapmap3.chr.8.pval.txt" 
     "/Users/christophergerman/.julia/dev/OrdinalGWAS/data/hapmap3.chr.9.pval.txt" 



In the rare situations where the multiple sets of Plink files lack the `fam` file or the corresponding bed and bim files have different filenames, users can explicitly supply bed filename, bim file name, and number of individuals. Replace Step 2 by 

Step 2': GWAS for each chromosome.


```julia
# this part can be submitted as separate jobs
for chr in 1:23
    bedfile = datadir * "hapmap3.chr." * string(chr) * ".bed"
    bimfile = datadir * "hapmap3.chr." * string(chr) * ".bim"
    pvalfile = datadir * "hapmap3.chr." * string(chr) * ".pval.txt"
    ordinalgwas(nm, bedfile, bimfile, 324; pvalfile=pvalfile)
end
```


```julia
# clean up
rm("ordinalgwas.null.txt", force=true)
isfile(datadir * "fittednullmodel.jld2") && rm(datadir * "fittednullmodel.jld2")
for chr in 1:23
    pvalfile = datadir * "hapmap3.chr." * string(chr) * ".pval.txt"
    rm(pvalfile, force=true)
end
for chr in 1:26
    plinkfile = datadir * "hapmap3.chr." * string(chr)
    rm(plinkfile * ".bed", force=true)
    rm(plinkfile * ".fam", force=true)
    rm(plinkfile * ".bim", force=true)
end
```

## Multiple Plink file sets on cluster

We provide two scripts that successfully run on UCLA's Hoffman2 cluster using Julia v1.0.1 and PBS job schedulaer (`qsub`).

* The first script [`cluster_preparedata.jl`](https://raw.githubusercontent.com/OpenMendel/OrdinalGWAS.jl/master/docs/cluster_preparedata.jl) creates a demo data set in current folder. Run
```julia
julia cluster_preparedata.jl
```
on head node.


```julia
run(`cat cluster_preparedata.jl`);
```

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


* The second script [`cluster_run.jl`](https://raw.githubusercontent.com/OpenMendel/OrdinalGWAS.jl/master/docs/cluster_run.jl) first fits the null model then submits a separate job for each chromosome. Run
```julia
julia cluster_run.jl
```
on head node.


```julia
run(`cat cluster_run.jl`);
```

    #!/usr/local/bin/julia
    #
    # This script demonstrates how to submit multiple OrdinalGWAS runs from multiple sets of
    # Plink files on UCLA Hoffman2 cluster. It assumes that a demo data is available by
    # running `julia cluster_preparedata.jl` at current folder.
    #
    
    using OrdinalGWAS, Serialization
    
    # Step 1: fit null model and save result to file `fittednullmodel.jls`
    nm = ordinalgwas(@formula(trait ~ sex), "covariate.txt", nothing)
    open("fittednullmodel.jls", "w") do io
        Serialization.serialize(io, nm)
    end
    
    # Step 2: GWAS for each chromosome
    for chr in 1:23
        println("submit job for chromosome=$chr")
        jcode = "using OrdinalGWAS, Serialization;
        nm = open(deserialize, \"fittednullmodel.jls\");
        bedfile = \"hapmap3.chr.\" * string($chr) * \".bed.gz\";
        bimfile = \"hapmap3.chr.\" * string($chr) * \".bim.gz\";
        pvalfile = \"hapmap3.chr.\" * string($chr) * \".pval.txt\";
        ordinalgwas(nm, bedfile, bimfile, 324; pvalfile=pvalfile);"
        # prepare sh file for qsub
        open("tmp.sh", "w") do io
            println(io, "#!/bin/bash")
            println(io, "#\$ -cwd")
            println(io, "# error = Merged with joblog")
            println(io, "#\$ -o joblog.\$JOB_ID")
            println(io, "#\$ -j y")
            println(io, "#\$ -l h_rt=0:30:00,h_data=2G") # request runtime and memory
            println(io, "#\$ -pe shared 2") # request # shared-memory nodes
            println(io, "# Email address to notify")
            println(io, "#\$ -M \$USER@mail")
            println(io, "# Notify when")
            println(io, "#\$ -m a")
            println(io)
            println(io, "# load the job environment:")
            println(io, ". /u/local/Modules/default/init/modules.sh")
            println(io, "module load julia/1.0.1") # available Julia version
            println(io)
            println(io, "# run julia code")
            println(io, "julia -e '$jcode' > output.\$JOB_ID 2>&1")
        end
        # submit job
        run(`qsub tmp.sh`)
    end



```julia

```
