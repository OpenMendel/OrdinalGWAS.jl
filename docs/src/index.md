
# PolrGWAS.jl

PolrGWAS.jl is a Julia package for performing genome-wide association studies (GWAS) for ordered categorical phenotypes using [proportional odds model](https://en.wikipedia.org/wiki/Ordered_logit) or [ordred Probit model](https://en.wikipedia.org/wiki/Ordered_probit). It is useful when the phenotype takes ordered discrete values, e.g., disease status (undiagnosed, pre-disease, mild, moderate, severe). The package name follows the function [`polr` in R package MASS](https://www.rdocumentation.org/packages/MASS/versions/7.3-3/topics/polr).

## Installation

This package requires Julia v0.7.0 or later and two other unregistered packages SnpArrays and PolrModels. The package has not yet been registered and must be installed using the repository location. Start julia and use the ] key to switch to the package manager REPL
```julia
(v0.7) pkg> add https://github.com/OpenMendel/SnpArrays.git#juliav0.7
(v0.7) pkg> add https://github.com/OpenMendel/PolrModels.git
(v0.7) pkg> add https://github.com/OpenMendel/PolrGWAS.git
```


```julia
# machine information for this tutorial
versioninfo()
```

    Julia Version 0.7.0
    Commit a4cb80f3ed (2018-08-08 06:46 UTC)
    Platform Info:
      OS: macOS (x86_64-apple-darwin14.5.0)
      CPU: Intel(R) Core(TM) i7-6920HQ CPU @ 2.90GHz
      WORD_SIZE: 64
      LIBM: libopenlibm
      LLVM: libLLVM-6.0.0 (ORCJIT, skylake)
    Environment:
      JULIA_EDITOR = code



```julia
# for use in this tutorial
using BenchmarkTools, CSV, PolrGWAS, SnpArrays
```

## Example data set

`data` folder of the package contains an example data set.


```julia
;ls -l ../data
```

    total 3664
    -rw-r--r--  1 huazhou  staff     6844 Nov 23 17:58 covariate.txt
    -rw-r--r--  1 huazhou  staff  1128171 Nov 23 17:58 hapmap3.bed
    -rw-r--r--  1 huazhou  staff   388672 Nov 23 17:58 hapmap3.bim
    -rw-r--r--  1 huazhou  staff     7136 Nov 23 17:58 hapmap3.fam
    -rw-r--r--  1 huazhou  staff   332960 Nov 23 17:58 hapmap3.map
    -rw-r--r--  1 huazhou  staff      773 Nov 23 17:58 simtrait.jl


`covariate.txt` is a comma separated value (CSV) file containing the sample information, covariates `sex`, and phenotype `trait`. `trait` is coded as integer values 1, 2, 3 or 4. It was simulated from the script `simtrait.jl`. 


```julia
;head -20 ../data/covariate.txt
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
    2426,NA19904,0,0,1,1
    2431,NA19917,0,0,2,1
    2436,NA19982,0,0,1,2
    2487,NA20340,0,0,1,4
    2427,NA19909,0,0,2,4
    2424,NA19834,0,0,1,4
    2480,NA20317,0,0,2,4
    2418,NA19818,0,0,1,1
    2490,NA20346,0,0,1,2
    2433,NA19921,0,0,2,4


`hapmap3` is a set of Plink files that contain the genotype information of samples. 

## Basic usage

The following command performs GWAS using the [proportional odds model](https://en.wikipedia.org/wiki/Ordered_logit).


```julia
polrgwas(@formula(trait ~ 0 + sex), "../data/covariate.txt", "../data/hapmap3")
```

For documentation of the `polrgwas` function, type `?polrgwas` in Julia REPL.
```@docs
polrgwas
```

### Formula for null model

The first argument specifies the null model without SNP effects, e.g., `@formula(trait ~ 0 + sex)`. Following R and Julia convention, `trait ~ 0 + sex` explicitly excludes intercept. It is important to **exclude** the intercept in `polr` model formulation because proportional odds model automatically incorporates intercepts for modeling purpose.

### Input files

`polrgwas` expects two input files: one for responses plus covariates (second argument), the other the Plink files for genotypes (third argument).

Covariates and phenotype are available in a csv file `covariate.txt`, which has one header line for variable names. Variable `trait` is the ordered categorical phenotypes coded as integers 1 to 4. We want to include variable `sex` as the covariate in GWAS.


```julia
;head -20 ../data/covariate.txt
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
    2426,NA19904,0,0,1,1
    2431,NA19917,0,0,2,1
    2436,NA19982,0,0,1,2
    2487,NA20340,0,0,1,4
    2427,NA19909,0,0,2,4
    2424,NA19834,0,0,1,4
    2480,NA20317,0,0,2,4
    2418,NA19818,0,0,1,1
    2490,NA20346,0,0,1,2
    2433,NA19921,0,0,2,4


Genotype data is available as binary Plink files.


```julia
;ls -l ../data/hapmap3.bed ../data/hapmap3.bim ../data/hapmap3.fam
```

    -rw-r--r--  1 huazhou  staff  1128171 Nov 23 17:58 ../data/hapmap3.bed
    -rw-r--r--  1 huazhou  staff   388672 Nov 23 17:58 ../data/hapmap3.bim
    -rw-r--r--  1 huazhou  staff     7136 Nov 23 17:58 ../data/hapmap3.fam


There are 324 samples at 13,928 SNPs.


```julia
size(SnpArray("../data/hapmap3.bed"))
```




    (324, 13928)



### Output files

`polrgwas` outputs two files: `polrgwas.nullmodel.txt` and `polrgwas.scoretest.txt`. The prefix `polrgwas` of output files can be changed by the `outfile` keyword, e.g.,
```julia
polrgwas(@formula(trait ~ 0 + sex), "../data/covariate.txt", "../data/hapmap3", 
    outfile="hapmap3")
```

* `polrgwas.nullmodel.txt` lists the estimated null model (without SNPs).  


```julia
;cat polrgwas.nullmodel.txt
```

    StatsModels.DataFrameRegressionModel{PolrModel{Int64,Float64,LogitLink},Array{Float64,2}}
    
    Formula: trait ~ +sex
    
    Coefficients:
          Estimate Std.Error  t value Pr(>|t|)
    θ1    -1.48564  0.358891 -4.13952    <1e-4
    θ2   -0.569479  0.341044 -1.66981   0.0959
    θ3    0.429815  0.339642  1.26549   0.2066
    β1    0.424656  0.213914  1.98517   0.0480


* `polrgwas.scoretest.txt` tallies the SNPs and their pvalues. 


```julia
;head polrgwas.scoretest.txt
```

    chr,pos,snpid,maf,pval
    1,554484,rs10458597,0.0,1.0
    1,758311,rs12562034,0.07763975155279501,0.004565312839540994
    1,967643,rs2710875,0.32407407407407407,3.1082838285548695e-5
    1,1168108,rs11260566,0.19158878504672894,1.2168672367668912e-5
    1,1375074,rs1312568,0.441358024691358,0.008206860046175221
    1,1588771,rs35154105,0.0,1.0
    1,1789051,rs16824508,0.00462962962962965,0.5111981332544
    1,1990452,rs2678939,0.4537037037037037,0.29972829571847825
    1,2194615,rs7553178,0.22685185185185186,0.1713331245805063


### Input non-genetic data as DataFrame

Internally `polrgwas` parses the covariate file as a DataFrame by `CSV.read(covfile)`. For covariate file of other format, users can parse first and then input a DataFrame to `polrgwas` directly.
```julia
polrgwas(@formula(trait ~ 0 + sex), df, "../data/hapmap3")
```
!!! note

    Users should always make sure that individuals in covariate file or DataFrame match those in Plink fam file. 

For example, following code checks that the first 2 columns of the `covariate.txt` file match the first 2 columns of the `hapmap3.fam` file exactly.


```julia
covdf = CSV.read("../data/covariate.txt")
plkfam = CSV.read("../data/hapmap3.fam", header=0, delim=' ')
all(covdf[1] .== plkfam[1]) && all(covdf[2] .== plkfam[2])
```




    true



### Timing

For this moderate-sized data set, `polrgwas` takes less than 0.2 second.


```julia
@btime(polrgwas(@formula(trait ~ 0 + sex), "../data/covariate.txt", "../data/hapmap3"))
```

      130.712 ms (639477 allocations: 29.13 MiB)



```julia
# clean up
rm("polrgwas.scoretest.txt")
rm("polrgwas.nullmodel.txt")
```

## Link functions

The `link` keyword argument of `polrgwas` can take value:  
- `LogitLink()`, proportional odds model (default),  
- `ProbitLink()`, ordred Probit model,  
- `CloglogLink()`, proportional hazards model), or 
- `CauchyLink()`.

For example, to perform GWAS using the ordred Probit model


```julia
polrgwas(@formula(trait ~ 0 + sex), "../data/covariate.txt", "../data/hapmap3", 
    link=ProbitLink(), outfile="opm")
```

The estimates in null model and p-values are slightly different from proportional odds moodel.


```julia
;cat opm.nullmodel.txt
```

    StatsModels.DataFrameRegressionModel{PolrModel{Int64,Float64,ProbitLink},Array{Float64,2}}
    
    Formula: trait ~ +sex
    
    Coefficients:
          Estimate Std.Error  t value Pr(>|t|)
    θ1   -0.866156  0.210677 -4.11129    <1e-4
    θ2   -0.359878  0.205817 -1.74854   0.0813
    θ3    0.247054  0.205382   1.2029   0.2299
    β1    0.251058  0.128225  1.95795   0.0511



```julia
;head opm.scoretest.txt
```

    chr,pos,snpid,maf,pval
    1,554484,rs10458597,0.0,1.0
    1,758311,rs12562034,0.07763975155279501,0.010076916742300138
    1,967643,rs2710875,0.32407407407407407,2.6272564941853933e-5
    1,1168108,rs11260566,0.19158878504672894,1.0897484851078458e-5
    1,1375074,rs1312568,0.441358024691358,0.005102883990438149
    1,1588771,rs35154105,0.0,1.0
    1,1789051,rs16824508,0.00462962962962965,0.48653776297859236
    1,1990452,rs2678939,0.4537037037037037,0.33231290090455434
    1,2194615,rs7553178,0.22685185185185186,0.25915513977197435



```julia
rm("opm.nullmodel.txt")
rm("opm.scoretest.txt")
```

## SNP models

Genotypes are translated into numeric values according to different genetic model, which is specified by the `snpmodel` keyword. Default is `ADDITIVE_MODEL`.

| Genotype | `SnpArray` | `ADDITIVE_MODEL` | `DOMINANT_MODEL` | `RECESSIVE_MODEL` |    
|:---:|:---:|:---:|:---:|:---:|  
| A1,A1 | 0x00 | 0 | 0 | 0 |  
| missing | 0x01 | NaN | NaN | NaN |
| A1,A2 | 0x02 | 1 | 1 | 0 |  
| A2,A2 | 0x03 | 2 | 1 | 1 |  

## SNP and/or sample masks

In practice, we often perform GWAS on selected SNPs and/or selected samples. They can be specified by the `colinds` and `rowinds` keywords of `polrgwas` function.

For example, to perform GWAS on SNPs with minor allele frequency (MAF) above 0.05


```julia
# create SNP mask
snpinds = maf(SnpArray("../data/hapmap3.bed")) .≥ 0.05
```




    13928-element BitArray{1}:
     false
      true
      true
      true
      true
     false
     false
      true
      true
      true
     false
      true
      true
         ⋮
      true
      true
      true
      true
      true
      true
      true
      true
      true
     false
     false
     false




```julia
# GWAS on selected SNPs
@time polrgwas(@formula(trait ~ 0 + sex), "../data/covariate.txt", "../data/hapmap3", 
    colinds = snpinds, outfile="commonvariant")
```

      0.219613 seconds (831.17 k allocations: 39.139 MiB, 3.44% gc time)



```julia
;cat commonvariant.nullmodel.txt
```

    StatsModels.DataFrameRegressionModel{PolrModel{Int64,Float64,LogitLink},Array{Float64,2}}
    
    Formula: trait ~ +sex
    
    Coefficients:
          Estimate Std.Error  t value Pr(>|t|)
    θ1    -1.48564  0.358891 -4.13952    <1e-4
    θ2   -0.569479  0.341044 -1.66981   0.0959
    θ3    0.429815  0.339642  1.26549   0.2066
    β1    0.424656  0.213914  1.98517   0.0480



```julia
;head -20 commonvariant.scoretest.txt
```

    chr,pos,snpid,maf,pval
    1,758311,rs12562034,0.07763975155279501,0.004565312839540994
    1,967643,rs2710875,0.32407407407407407,3.1082838285548695e-5
    1,1168108,rs11260566,0.19158878504672894,1.2168672367668912e-5
    1,1375074,rs1312568,0.441358024691358,0.008206860046175221
    1,1990452,rs2678939,0.4537037037037037,0.29972829571847825
    1,2194615,rs7553178,0.22685185185185186,0.1713331245805063
    1,2396747,rs13376356,0.1448598130841121,0.5320416198875456
    1,2823603,rs1563468,0.4830246913580247,0.225191391783573
    1,3025087,rs6690373,0.2538699690402477,0.7018469417717486
    1,3431124,rs12093117,0.1099071207430341,0.42779374684900734
    1,3633945,rs10910017,0.22187500000000004,0.913128778096905
    1,4096895,rs6702633,0.4752321981424149,0.006516308271075796
    1,4297388,rs684965,0.3055555555555556,0.09519531577509736
    1,4498133,rs11809295,0.0993788819875776,0.08324350470610889
    1,4698713,rs578528,0.32407407407407407,0.06923073232117855
    1,4899946,rs4654471,0.3580246913580247,0.22453129666689242
    1,5100369,rs6681148,0.13157894736842102,0.15566709962762512
    1,5302730,rs10799197,0.4287925696594427,0.6690547960840125
    1,5502779,rs10796400,0.2314814814814815,0.24152455743940068



```julia
# extra headline in commonvariant.scoretest.txt
countlines("commonvariant.scoretest.txt"), count(snpinds)
```




    (12086, 12085)




```julia
# clean up
rm("commonvariant.scoretest.txt")
rm("commonvariant.nullmodel.txt")
```

User should be particularly careful when using the `rowinds` keyword. Selected rows in SnpArray should exactly match the samples in the null model. Otherwise the results are meaningless.

## Likelihood ratio test (LRT)

By default, `polrgwas` calculates p-value for each SNP using score test. Score test is fast because it doesn't require fitting alternative model for each SNP. User can request likelihood ratio test (LRT) using keyword `test=:LRT`. LRT is much slower but may be more powerful than score test.


```julia
@time polrgwas(@formula(trait ~ 0 + sex), "../data/covariate.txt", "../data/hapmap3", 
    test=:LRT, outfile="lrt")
```

     20.709661 seconds (8.81 M allocations: 2.065 GiB, 1.80% gc time)


Test result is output to `outfile.lrttest.txt` file


```julia
;head -20 lrt.lrttest.txt
```

    chr,pos,snpid,maf,effect,pval
    1,554484,rs10458597,0.0,0.0,1.0
    1,758311,rs12562034,0.07763975155279501,-1.0057833719544331,0.0019185836579804134
    1,967643,rs2710875,0.32407407407407407,-0.6488560566295055,1.805050556976241e-5
    1,1168108,rs11260566,0.19158878504672894,-0.9157225669357879,5.873384712685568e-6
    1,1375074,rs1312568,0.441358024691358,-0.33181366525772593,0.008081022577832324
    1,1588771,rs35154105,0.0,0.0,1.0
    1,1789051,rs16824508,0.00462962962962965,-0.7338026388701573,0.5169027130129711
    1,1990452,rs2678939,0.4537037037037037,-0.13586499231819726,0.29946402200912603
    1,2194615,rs7553178,0.22685185185185186,-0.2512075640440123,0.16151069094439868
    1,2396747,rs13376356,0.1448598130841121,0.12946142026273783,0.5387338201469207
    1,2623158,rs28753913,0.0,0.0,1.0
    1,2823603,rs1563468,0.4830246913580247,0.15515329587697405,0.23123002081577315
    1,3025087,rs6690373,0.2538699690402477,-0.05966638389967704,0.6995722170701131
    1,3225416,rs12043519,0.029320987654321007,1.1761887120778427,0.002016744167744886
    1,3431124,rs12093117,0.1099071207430341,0.1824233299545801,0.4305201231293997
    1,3633945,rs10910017,0.22187500000000004,-0.01793569262704999,0.9142024828415289
    1,3895935,rs34770924,0.024691358024691357,0.0009448575482649847,0.9980482909192042
    1,4096895,rs6702633,0.4752321981424149,0.42308741025632207,0.0063052493446845445
    1,4297388,rs684965,0.3055555555555556,-0.27091872232225434,0.09226258176075279


Note the extra `effect` column, which is the effect size (regression coefficient) for each SNP. 


```julia
# clean up
rm("lrt.lrttest.txt")
rm("lrt.nullmodel.txt")
```

In this example, GWAS by score test takes less than 0.2 second, while GWAS by LRT takes about 20 seconds. About 100 fold difference in run time. 

## Score test for screening, LRT for power 

For large data sets, a practical solution is to perform score test first, then re-do LRT for the most promising SNPs according to score test p-values.

**Step 1**: Perform score test GWAS, results in `hapmap.scoretest.txt`.


```julia
@time polrgwas(@formula(trait ~ 0 + sex), "../data/covariate.txt", "../data/hapmap3", 
    test=:score, outfile="hapmap", verbose=false)
```

      0.249964 seconds (718.42 k allocations: 33.304 MiB, 12.06% gc time)



```julia
;head -20 hapmap.scoretest.txt
```

    chr,pos,snpid,maf,pval
    1,554484,rs10458597,0.0,1.0
    1,758311,rs12562034,0.07763975155279501,0.004565312839540994
    1,967643,rs2710875,0.32407407407407407,3.1082838285548695e-5
    1,1168108,rs11260566,0.19158878504672894,1.2168672367668912e-5
    1,1375074,rs1312568,0.441358024691358,0.008206860046175221
    1,1588771,rs35154105,0.0,1.0
    1,1789051,rs16824508,0.00462962962962965,0.5111981332544
    1,1990452,rs2678939,0.4537037037037037,0.29972829571847825
    1,2194615,rs7553178,0.22685185185185186,0.1713331245805063
    1,2396747,rs13376356,0.1448598130841121,0.5320416198875456
    1,2623158,rs28753913,0.0,1.0
    1,2823603,rs1563468,0.4830246913580247,0.225191391783573
    1,3025087,rs6690373,0.2538699690402477,0.7018469417717486
    1,3225416,rs12043519,0.029320987654321007,0.0010797837613683529
    1,3431124,rs12093117,0.1099071207430341,0.42779374684900734
    1,3633945,rs10910017,0.22187500000000004,0.913128778096905
    1,3895935,rs34770924,0.024691358024691357,0.9990210608648795
    1,4096895,rs6702633,0.4752321981424149,0.006516308271075796
    1,4297388,rs684965,0.3055555555555556,0.09519531577509736


**Step 2**: Sort score test p-values and find top 10 SNPs.


```julia
scorepvals = CSV.read("hapmap.scoretest.txt")[5]
```




    13928-element Array{Union{Missing, Float64},1}:
     1.0                  
     0.004565312839540994 
     3.1082838285548695e-5
     1.2168672367668912e-5
     0.008206860046175221 
     1.0                  
     0.5111981332544      
     0.29972829571847825  
     0.1713331245805063   
     0.5320416198875456   
     1.0                  
     0.225191391783573    
     0.7018469417717486   
     ⋮                    
     0.3659920717810955   
     0.5612543378640146   
     0.11646522388490779  
     0.5216852471858483   
     0.44711860378438606  
     0.16565527326637602  
     0.65590930553447     
     0.89481895007441     
     0.6559093055344556   
     0.3462815220311233   
     0.2316733923395225   
     0.2651082379376032   




```julia
tophits = sortperm(scorepvals)[1:10]
```




    10-element Array{Int64,1}:
      6063
      5071
      2458
         4
      3291
     13737
      5293
         3
      6256
      4183




```julia
scorepvals[tophits]
```




    10-element Array{Union{Missing, Float64},1}:
     1.3080149099181335e-6
     6.536722765052079e-6 
     9.664742185669054e-6 
     1.2168672367668912e-5
     1.802746001833127e-5 
     2.0989542284213636e-5
     2.6844521269963608e-5
     3.1082838285548695e-5
     4.1010912875160476e-5
     4.2966265138454806e-5



**Step 3**: Re-do LRT on top hits.


```julia
@time polrgwas(@formula(trait ~ 0 + sex), "../data/covariate.txt", "../data/hapmap3", 
    colinds=tophits, test=:LRT, outfile="hapmap")
```

      0.162109 seconds (411.63 k allocations: 21.434 MiB, 4.57% gc time)



```julia
;head -20 hapmap.lrttest.txt
```

    chr,pos,snpid,maf,effect,pval
    1,967643,rs2710875,0.32407407407407407,-0.6488560566295055,1.805050556976241e-5
    1,1168108,rs11260566,0.19158878504672894,-0.9157225669357879,5.873384712685568e-6
    3,36821790,rs4678553,0.23456790123456794,0.7424952268973518,1.1303825016262592e-5
    4,11017683,rs16881446,0.27554179566563464,-0.7870581482955515,1.1105427468799613e-5
    5,3739190,rs12521166,0.0679012345679012,1.1468852997925316,4.781288229657399e-5
    6,7574576,rs1885466,0.17746913580246915,0.8750621092263019,7.272346896740631e-6
    6,52474721,rs2073183,0.1826625386996904,0.7790794914858663,5.069394513906121e-5
    7,41152376,rs28880,0.3379629629629629,-0.814633902445351,9.180126530294943e-7
    7,84223996,rs4128623,0.07870370370370372,1.0022229316338573,6.587895464657512e-5
    23,121048059,rs1937165,0.4380804953560371,0.5392313636256612,1.9754643855522616e-5



```julia
# clean up
rm("hapmap.nullmodel.txt")
rm("hapmap.lrttest.txt")
rm("hapmap.scoretest.txt")
```

## GxE or other interactions

In many applications, we want to test SNP effect and/or its interaction with other terms. `testformula` keyword specifies the test unit **besides** the covariates in `nullformula`. 

In following example, keyword `testformula=@formula(trait ~ 0 + snp + snp & sex)` instructs `polrgwas` to test joint effect of `snp` and `snp & sex` interaction.


```julia
polrgwas(@formula(trait ~ 0 + sex), "../data/covariate.txt", "../data/hapmap3", 
    outfile="GxE", testformula=@formula(trait ~ 0 + snp + snp & sex))
```


```julia
;head -20 GxE.scoretest.txt
```

    chr,pos,snpid,maf,pval
    1,554484,rs10458597,0.0,1.0
    1,758311,rs12562034,0.07763975155279501,0.017446010412254197
    1,967643,rs2710875,0.32407407407407407,0.0001667073239489097
    1,1168108,rs11260566,0.19158878504672894,4.763762457893366e-5
    1,1375074,rs1312568,0.441358024691358,0.029138471242993652
    1,1588771,rs35154105,0.0,1.0
    1,1789051,rs16824508,0.00462962962962965,0.2964363114944328
    1,1990452,rs2678939,0.4537037037037037,0.37924580479348785
    1,2194615,rs7553178,0.22685185185185186,0.325582269932396
    1,2396747,rs13376356,0.1448598130841121,0.81664019848141
    1,2623158,rs28753913,0.0,1.0
    1,2823603,rs1563468,0.4830246913580247,0.394387712136073
    1,3025087,rs6690373,0.2538699690402477,0.929341576466736
    1,3225416,rs12043519,0.029320987654321007,0.004580867523699672
    1,3431124,rs12093117,0.1099071207430341,0.521577659391361
    1,3633945,rs10910017,0.22187500000000004,0.9879831700595433
    1,3895935,rs34770924,0.024691358024691357,0.39083539216130736
    1,4096895,rs6702633,0.4752321981424149,0.012628962362857694
    1,4297388,rs684965,0.3055555555555556,0.15091465952476424



```julia
# clean up
rm("GxE.nullmodel.txt")
rm("GxE.scoretest.txt")
```

To test the linear and quadratic SNP effects jointly


```julia
polrgwas(@formula(trait ~ 0 + sex), "../data/covariate.txt", "../data/hapmap3", 
    outfile="quadratic", testformula=@formula(trait ~ 0 + snp + snp & snp))
```


```julia
;head -20 quadratic.scoretest.txt
```

    chr,pos,snpid,maf,pval
    1,554484,rs10458597,0.0,1.0
    1,758311,rs12562034,0.07763975155279501,1.0
    1,967643,rs2710875,0.32407407407407407,1.0
    1,1168108,rs11260566,0.19158878504672894,1.0
    1,1375074,rs1312568,0.441358024691358,1.0
    1,1588771,rs35154105,0.0,1.0
    1,1789051,rs16824508,0.00462962962962965,1.0
    1,1990452,rs2678939,0.4537037037037037,1.0
    1,2194615,rs7553178,0.22685185185185186,1.0
    1,2396747,rs13376356,0.1448598130841121,1.0
    1,2623158,rs28753913,0.0,1.0
    1,2823603,rs1563468,0.4830246913580247,1.0
    1,3025087,rs6690373,0.2538699690402477,1.0
    1,3225416,rs12043519,0.029320987654321007,1.0
    1,3431124,rs12093117,0.1099071207430341,1.0
    1,3633945,rs10910017,0.22187500000000004,1.0
    1,3895935,rs34770924,0.024691358024691357,0.9999992473357002
    1,4096895,rs6702633,0.4752321981424149,1.0
    1,4297388,rs684965,0.3055555555555556,1.0



```julia
# clean up
rm("quadratic.nullmodel.txt")
rm("quadratic.scoretest.txt")
```

## Plots

To plot the GWAS results, use the MendelPlots (TODO) package.

## Docker

For ease of using PolrGWAS, we provide a Dockerfile so users don't need to install Julia and required packages. Only Docker app needs to be installed in order to run analysis. Following is tested on Docker 2.0.0.0-mac78.

**Step 1**: Create a Dockerfile with content shown below. You can copy from [here](https://raw.githubusercontent.com/OpenMendel/PolrGWAS.jl/master/docker/Dockerfile).


```julia
;cat ../docker/Dockerfile
```

    FROM centos:7
    
    WORKDIR /root
    
    RUN yum update -y && yum install -y epel-release && yum clean all
    
    RUN yum update -y && yum install -y \
        cmake \
        curl-devel \
        expat-devel \
        gcc \
        gcc-c++ \
        gcc-gfortran \
        gettext-devel \
        make \
        openssl \
        openssl098e \
        openssl-devel \
        patch \
        svn \
        wget \
        which \
        && yum clean all    
    
    ENV PATH /usr/local/sbin:/usr/local/bin:$PATH
    
    ENV LD_LIBRARY_PATH /usr/local/lib:/usr/local/lib64
    
    # GIT - https://git-scm.com/
    # http://tecadmin.net/install-git-2-0-on-centos-rhel-fedora/#
    ENV GIT_VER 2.19.1
    
    RUN wget https://www.kernel.org/pub/software/scm/git/git-$GIT_VER.tar.gz \
        && tar xf git-$GIT_VER.tar.gz && cd git-$GIT_VER \
        && make -j"$(nproc --all)" prefix=/usr/local all \
        && make prefix=/usr/local -j"$(nproc --all)" install \
        && cd .. && rm -f git-$GIT_VER.tar.gz && rm -rf git-$GIT_VER
    
    # Makes git use https by default
    RUN git config --global url."https://".insteadOf git://
    
    # Julia
    ENV JULIA_VER_MAJ 1.0
    ENV JULIA_VER_MIN .1
    ENV JULIA_VER $JULIA_VER_MAJ$JULIA_VER_MIN
    
    RUN wget https://julialang-s3.julialang.org/bin/linux/x64/$JULIA_VER_MAJ/julia-$JULIA_VER-linux-x86_64.tar.gz \
        && mkdir /usr/local/julia \
        && tar xf julia-$JULIA_VER-linux-x86_64.tar.gz --directory /usr/local/julia --strip-components=1 \
        && ln -s /usr/local/julia/bin/julia /usr/local/bin/julia \
        && rm -f julia-$JULIA_VER-linux-x86_64.tar.gz
    
    ENV JULIA_PKGDIR /usr/local/julia/share/julia/site
    
    RUN julia -e 'using Pkg; \
        Pkg.add([ \
        PackageSpec(url="https://github.com/OpenMendel/SnpArrays.jl.git", rev="juliav0.7"), \
        PackageSpec(url="https://github.com/OpenMendel/PolrModels.jl.git", rev="master"), \
        PackageSpec(url="https://github.com/OpenMendel/PolrGWAS.jl.git", rev="master") \
        ]); \
        Pkg.test("PolrGWAS");'


**Step 2**: Build a docker image `polrgwas-app`, assuming that the Dockerfile is located at `../docker` folder. Building the image for the first time can take up to 10 minutes; but it only needs to be done once.


```julia
;docker build -t polrgwas-app ../docker/
```

    Sending build context to Docker daemon  3.584kB
    Step 1/15 : FROM centos:7
     ---> 75835a67d134
    Step 2/15 : WORKDIR /root
     ---> Using cache
     ---> f052f835f42d
    Step 3/15 : RUN yum update -y && yum install -y epel-release && yum clean all
     ---> Using cache
     ---> a75b6baa5581
    Step 4/15 : RUN yum update -y && yum install -y     cmake     curl-devel     expat-devel     gcc     gcc-c++     gcc-gfortran     gettext-devel     make     openssl     openssl098e     openssl-devel     patch     svn     wget     which     && yum clean all
     ---> Using cache
     ---> dcf0f6ba6334
    Step 5/15 : ENV PATH /usr/local/sbin:/usr/local/bin:$PATH
     ---> Using cache
     ---> 7861f46c17eb
    Step 6/15 : ENV LD_LIBRARY_PATH /usr/local/lib:/usr/local/lib64
     ---> Using cache
     ---> 4d3222181efc
    Step 7/15 : ENV GIT_VER 2.19.1
     ---> Using cache
     ---> b23c5cb67786
    Step 8/15 : RUN wget https://www.kernel.org/pub/software/scm/git/git-$GIT_VER.tar.gz     && tar xf git-$GIT_VER.tar.gz && cd git-$GIT_VER     && make -j"$(nproc --all)" prefix=/usr/local all     && make prefix=/usr/local -j"$(nproc --all)" install     && cd .. && rm -f git-$GIT_VER.tar.gz && rm -rf git-$GIT_VER
     ---> Using cache
     ---> 6beaafd797ec
    Step 9/15 : RUN git config --global url."https://".insteadOf git://
     ---> Using cache
     ---> f1aec5dd453f
    Step 10/15 : ENV JULIA_VER_MAJ 1.0
     ---> Using cache
     ---> 42eb17d93684
    Step 11/15 : ENV JULIA_VER_MIN .1
     ---> Using cache
     ---> 4289e3f44a22
    Step 12/15 : ENV JULIA_VER $JULIA_VER_MAJ$JULIA_VER_MIN
     ---> Using cache
     ---> 2e9b9213c88c
    Step 13/15 : RUN wget https://julialang-s3.julialang.org/bin/linux/x64/$JULIA_VER_MAJ/julia-$JULIA_VER-linux-x86_64.tar.gz     && mkdir /usr/local/julia     && tar xf julia-$JULIA_VER-linux-x86_64.tar.gz --directory /usr/local/julia --strip-components=1     && ln -s /usr/local/julia/bin/julia /usr/local/bin/julia     && rm -f julia-$JULIA_VER-linux-x86_64.tar.gz
     ---> Using cache
     ---> 2510e817738a
    Step 14/15 : ENV JULIA_PKGDIR /usr/local/julia/share/julia/site
     ---> Using cache
     ---> 3692e40f50fa
    Step 15/15 : RUN julia -e 'using Pkg;     Pkg.add([     PackageSpec(url="https://github.com/OpenMendel/SnpArrays.jl.git", rev="juliav0.7"),     PackageSpec(url="https://github.com/OpenMendel/PolrModels.jl.git", rev="master"),     PackageSpec(url="https://github.com/OpenMendel/PolrGWAS.jl.git", rev="master")     ]);     Pkg.test("PolrGWAS");'
     ---> Using cache
     ---> af7e8765610c
    Successfully built af7e8765610c
    Successfully tagged polrgwas-app:latest


**Step 3**: Suppose data files are located at `/Users/huazhou/.julia/dev/PolrGWAS/data` folder, run analysis by


```julia
;docker run -v /Users/huazhou/.julia/dev/PolrGWAS/data:/data -t polrgwas-app julia -e 'using PolrGWAS; polrgwas(@formula(trait ~ 0 + sex), "/data/covariate.txt", "/data/hapmap3", outfile="/data/polrgwas");'
```

Here  
- `-t polrgwas-app` creates a container using the `polrgwas-app` image build in step 2.  
- `-v /Users/huazhou/.julia/dev/PolrGWAS/data:/data` maps the `/Users/huazhou/.julia/dev/PolrGWAS/data` folder on host machine to the `/data` folder within the container. 
- `julia -e 'using PolrGWAS; polrgwas(@formula(trait ~ 0 + sex), "/data/covariate.txt", "/data/hapmap3", outfile="/data/polrgwas");` calls Julia and runs `polrgwas` function. 

The output files are at `/Users/huazhou/.julia/dev/PolrGWAS/data`. 


```julia
;ls -l /Users/huazhou/.julia/dev/PolrGWAS/data/
```

    total 5320
    -rw-r--r--  1 huazhou  staff     6844 Nov 23 17:58 covariate.txt
    -rw-r--r--  1 huazhou  staff  1128171 Nov 23 17:58 hapmap3.bed
    -rw-r--r--  1 huazhou  staff   388672 Nov 23 17:58 hapmap3.bim
    -rw-r--r--  1 huazhou  staff     7136 Nov 23 17:58 hapmap3.fam
    -rw-r--r--  1 huazhou  staff   332960 Nov 23 17:58 hapmap3.map
    -rw-r--r--  1 huazhou  staff      347 Dec  1 10:12 polrgwas.nullmodel.txt
    -rw-r--r--  1 huazhou  staff   842401 Dec  1 10:12 polrgwas.scoretest.txt
    -rw-r--r--  1 huazhou  staff      773 Nov 23 17:58 simtrait.jl



```julia
# clean up
rm("/Users/huazhou/.julia/dev/PolrGWAS/data/polrgwas.nullmodel.txt")
rm("/Users/huazhou/.julia/dev/PolrGWAS/data/polrgwas.scoretest.txt")
```
