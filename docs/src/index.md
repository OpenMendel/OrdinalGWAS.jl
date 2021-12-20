# OrdinalGWAS.jl

OrdinalGWAS.jl is a Julia package for performing genome-wide association studies (GWAS) for ordered categorical phenotypes using [proportional odds model](https://en.wikipedia.org/wiki/Ordered_logit) or [ordred Probit model](https://en.wikipedia.org/wiki/Ordered_probit). It is useful when the phenotype takes ordered discrete values, e.g., disease status (undiagnosed, pre-disease, mild, moderate, severe). 

The methods and applications of this software package are detailed in the following publication:

*German CA, Sinsheimer JS, Klimentidis YC, Zhou H, Zhou JJ. (2020) Ordered multinomial regression for genetic association analysis of ordinal phenotypes at Biobank scale. Genetic Epidemiology. 44:248-260. [https://doi.org/10.1002/gepi.22276](https://doi.org/10.1002/gepi.22276)*

OrdinalGWAS.jl currently supports [PLINK](https://zzz.bwh.harvard.edu/plink/), [VCF](https://en.wikipedia.org/wiki/Variant_Call_Format) (both dosage and genotype data) file formats, and [BGEN](https://www.well.ox.ac.uk/~gav/bgen_format/) file formats. We plan to add [PGEN](https://www.cog-genomics.org/plink/2.0/formats#pgen) support in the future. 

## Installation

This package requires Julia v1.6 or later. The package has been registered and can be installed using Julia package manager. Start julia and type:
```julia
using Pkg
Pkg.add("OrdinalGWAS")
```


```julia
# machine information for this tutorial
versioninfo()
```

    Julia Version 1.6.2
    Commit 1b93d53fc4 (2021-07-14 15:36 UTC)
    Platform Info:
      OS: macOS (x86_64-apple-darwin18.7.0)
      CPU: Intel(R) Core(TM) i7-7820HQ CPU @ 2.90GHz
      WORD_SIZE: 64
      LIBM: libopenlibm
      LLVM: libLLVM-11.0.1 (ORCJIT, skylake)



```julia
# for use in this tutorial
using BenchmarkTools, CSV, Glob, SnpArrays, OrdinalGWAS, DataFrames
```

## Example data sets

The `data` folder of the package contains the example data sets for use with PLINK and VCF Files. In general, the user can locate this folder by command:


```julia
using OrdinalGWAS
const datadir = normpath(joinpath(dirname(pathof(OrdinalGWAS)), "../data/"))
```




    "/Users/xyz/.julia/dev/OrdinalGWAS/data/"




```julia
# content of the data folder
readdir(glob"*.*", datadir)
```




    16-element Vector{String}:
     "/Users/xyz/.julia/dev/OrdinalGWAS/data/bgen_ex.csv"
     "/Users/xyz/.julia/dev/OrdinalGWAS/data/bgen_snpsetfile.txt"
     "/Users/xyz/.julia/dev/OrdinalGWAS/data/bgen_test.bgen"
     "/Users/xyz/.julia/dev/OrdinalGWAS/data/bgen_test.bgen.bgi"
     "/Users/xyz/.julia/dev/OrdinalGWAS/data/covariate.txt"
     "/Users/xyz/.julia/dev/OrdinalGWAS/data/hapmap3.bed"
     "/Users/xyz/.julia/dev/OrdinalGWAS/data/hapmap3.bim"
     "/Users/xyz/.julia/dev/OrdinalGWAS/data/hapmap3.fam"
     "/Users/xyz/.julia/dev/OrdinalGWAS/data/hapmap3.map"
     "/Users/xyz/.julia/dev/OrdinalGWAS/data/hapmap_snpsetfile.txt"
     "/Users/xyz/.julia/dev/OrdinalGWAS/data/simtrait.jl"
     "/Users/xyz/.julia/dev/OrdinalGWAS/data/simtrait_bgen.jl"
     "/Users/xyz/.julia/dev/OrdinalGWAS/data/simtrait_vcf.jl"
     "/Users/xyz/.julia/dev/OrdinalGWAS/data/snpsetfile_vcf.txt"
     "/Users/xyz/.julia/dev/OrdinalGWAS/data/vcf_example.csv"
     "/Users/xyz/.julia/dev/OrdinalGWAS/data/vcf_test.vcf.gz"



The `hapmap3` files `covariate.txt` file correspond to data examples using PLINK formatted files (.bed, .bim, .fam). 

The `vcf_test.vcf.gz` and `vcf_example.csv` files are for an example analysis using VCF formatted files. 


## Basic usage

The following command performs GWAS using the [proportional odds model](https://en.wikipedia.org/wiki/Ordered_logit) for the hapmap3 PLINK files. The output is the fitted null model.

The default type of GWAS performed is a single-snp significance genome-wide scan, this can be changed by the keyword `analysistype` (default is "singlesnp"). Other types of analyses are gone over later. 


```julia
ordinalgwas(@formula(trait ~ sex), datadir * "covariate.txt", datadir * "hapmap3")
```




    StatsModels.TableRegressionModel{OrdinalMultinomialModel{Int64, Float64, LogitLink}, Matrix{Float64}}
    
    trait ~ sex
    
    Coefficients:
    ──────────────────────────────────────────────────────
                   Estimate  Std.Error   t value  Pr(>|t|)
    ──────────────────────────────────────────────────────
    intercept1|2  -1.48564    0.358891  -4.13952    <1e-04
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

`ordinalgwas` expects two input files: one for responses plus covariates (second argument), the other the genetic file(s) for dosages/genotypes (third argument).

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


#### Genetic file(s)

OrdinalGWAS supports PLINK, VCF, and BGEN Files.

- Genotype data is available as binary PLINK files.

- OrdinalGWAS can use dosage or genotype data from VCF Files. 

- OrdinalGWAS can use dosage data from BGEN files.   

!!! note

    By default, OrdinalGWAS assumes a set of PLINK files will be used. When using a VCF File, VCF file and type of data (dosage, genotype) must be specified by the `geneticformat` and `vcftype` options (as shown later). Similarly, BGEN must be specified as the `geneticformat` if using a BGEN file. 


```julia
readdir(glob"hapmap3.*", datadir)
```




    4-element Vector{String}:
     "/Users/xyz/.julia/dev/OrdinalGWAS/data/hapmap3.bed"
     "/Users/xyz/.julia/dev/OrdinalGWAS/data/hapmap3.bim"
     "/Users/xyz/.julia/dev/OrdinalGWAS/data/hapmap3.fam"
     "/Users/xyz/.julia/dev/OrdinalGWAS/data/hapmap3.map"



In this example, there are 324 samples at 13,928 SNPs.


```julia
size(SnpArray(datadir * "hapmap3.bed"))
```




    (324, 13928)



Compressed PLINK and VCF files are supported. For example, if Plink files are `hapmap3.bed.gz`, `hapmap3.bim.gz` and `hapmap3.fam.gz`, the same command
```julia
ordinalgwas(@formula(trait ~ sex), datadir * "covariate.txt", datadir * "hapmap3")
```
still works. Check all supported compression format by


```julia
SnpArrays.ALLOWED_FORMAT
```




    6-element Vector{String}:
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

    StatsModels.TableRegressionModel{OrdinalMultinomialModel{Int64, Float64, LogitLink}, Matrix{Float64}}
    
    trait ~ sex
    
    Coefficients:
    ──────────────────────────────────────────────────────
                   Estimate  Std.Error   t value  Pr(>|t|)
    ──────────────────────────────────────────────────────
    intercept1|2  -1.48564    0.358891  -4.13952    <1e-04
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
    1,758311,rs12562034,0.07763975155279501,0.4098763332666681,0.004565312839535881
    1,967643,rs2710875,0.32407407407407407,4.076249100705747e-7,3.108283828553597e-5
    1,1168108,rs11260566,0.19158878504672894,0.1285682279446898,1.2168672367654906e-5
    1,1375074,rs1312568,0.441358024691358,2.5376019650614977e-19,0.008206860046174836
    1,1588771,rs35154105,0.0,1.0,1.0
    1,1789051,rs16824508,0.00462962962962965,0.9332783156468178,0.5111981332534471
    1,1990452,rs2678939,0.4537037037037037,5.07695957708431e-11,0.29972829571847304
    1,2194615,rs7553178,0.22685185185185186,0.17056143157457776,0.17133312458049288


Output file names can be changed by the `nullfile` and `pvalfile` keywords respectively. For example, 
```julia
ordinalgwas(@formula(trait ~ sex), datadir * "covariate.txt", datadir * "hapmap3",
    pvalfile="ordinalgwas.pval.txt.gz")
```
will output the p-value file in compressed gz format.

### Subsamples

Use the keyword `covrowinds` to specify selected samples in the covarite file. Use the keyword `geneticrowinds` to specify selected samples in the Plink bed file or VCF File. For example, to use the first 300 samples in both covariate and bed file:
```julia
ordinalgwas(@formula(trait ~ sex), covfile, geneticfile, 
    covrowinds=1:300, geneticrowinds=1:300)
```
!!! note

    Users should always make sure that the selected samples in covariate file match exactly those in bed file. 

### Input non-genetic data as DataFrame

Internally `ordinalgwas` parses the covariate file as a DataFrame by `CSV.read(covfile)`. For covariate file of other formats, users can parse it as a DataFrame and then input the DataFrame to `ordinalgwas` directly.
```julia
ordinalgwas(@formula(trait ~ sex), df, geneticfile)
```
!!! note

    Users should always make sure that individuals in covariate file or DataFrame match those in Plink fam file/VCF File/BGEN file. 

For example, following code checks that the first 2 columns of the `covariate.txt` file match the first 2 columns of the `hapmap3.fam` file exactly.


```julia
covdf = CSV.read(datadir * "covariate.txt", DataFrame)
plkfam = CSV.read(datadir * "hapmap3.fam", header=0, delim=' ',  DataFrame)
all(covdf[!, 1] .== plkfam[!, 1]) && all(covdf[!, 2] .== plkfam[!, 2])
```




    true



### Timing

For this moderate-sized data set, `ordinalgwas` takes around 0.2 seconds.


```julia
@btime(ordinalgwas(@formula(trait ~ sex), datadir * "covariate.txt", datadir * "hapmap3"));
```

      201.823 ms (500794 allocations: 49.17 MiB)



```julia
# clean up
rm("ordinalgwas.null.txt", force=true)
rm("ordinalgwas.pval.txt", force=true)
```

## VCF Formatted Files

By default, OrdinalGWAS.jl will assume you are using PLINK files. It also supports VCF (and BGEN) Files. `vcf.gz` files are supported as well. To use vcf files in any of the analysis options detailed in this documentation, you simply need to add two keyword options to the `ordinalgwas` function:
* `geneticformat`: Choices are "VCF" or "PLINK". If you are using a VCF file, use `geneticformat = "VCF"`.
* `vcftype`: Choices are :GT (for genotypes) or :DS (for dosages). This tells OrdinalGWAS which type of data to extract from the VCF file.

Using a VCF File does not output minor allele frequency or hardy weinberg equillibrium p-values for each SNP tested since they may be dosages and MAF/HWE calculations via VCFTools does not currently support dosage data. 

The following shows how to run an analysis with a VCF file using the dosage information. 


```julia
ordinalgwas(@formula(y ~ sex), datadir * "vcf_example.csv", 
    datadir * "vcf_test"; geneticformat = "VCF", vcftype = :DS)
```




    StatsModels.TableRegressionModel{OrdinalMultinomialModel{Int64, Float64, LogitLink}, Matrix{Float64}}
    
    y ~ sex
    
    Coefficients:
    ───────────────────────────────────────────────────────
                   Estimate  Std.Error    t value  Pr(>|t|)
    ───────────────────────────────────────────────────────
    intercept1|2  -1.10106    0.205365  -5.36147     <1e-06
    intercept2|3   0.370894   0.188822   1.96425     0.0510
    intercept3|4   1.74736    0.232756   7.50726     <1e-11
    sex            0.22796    0.262237   0.869289    0.3858
    ───────────────────────────────────────────────────────




```julia
run(`head ordinalgwas.pval.txt`);
```

    chr,pos,snpid,pval
    22,20000086,rs138720731,0.6435068072069544
    22,20000146,rs73387790,1.0
    22,20000199,rs183293480,0.9378258278500581
    22,20000291,rs185807825,0.21907288710091172
    22,20000428,rs55902548,0.002790402446884929
    22,20000683,rs142720028,1.0
    22,20000771,rs114690707,0.23034910158749672
    22,20000793,rs189842693,0.15612228776953083
    22,20000810,rs147349046,0.1741386270145462


## BGEN Formatted Files

By default, OrdinalGWAS.jl will assume you are using PLINK files. It also supports BGEN Files. To use BGEN files in any of the analysis options detailed in this documentation, you simply need to add the following keyword option to the `ordinalgwas` function:
* `geneticformat`: Choices are "VCF" or "PLINK" or "BGEN". If you are using a BGEN file, use `geneticformat = "BGEN"`.

Some features, such as SNP-set analyses, are only available currently when there's an indexing `.bgi` file available. 

!!! note

    BGEN files can contain an optional index file (`.bgi` file) that allows the variants to be specified in order of position. OrdinalGWAS will automatically look for a file in the same directory as the `BGENFILENAME` with the name `BGENFILENAME.bgi`. The BGEN file is read either in the `.bgi` file order if `BGENFILENAME.bgi` is supplied in the same directory as the BGEN file, otherwise it will use the order in the BGEN file. This is important in analyses specifying `snpinds` as well as annotation groupings. You must make sure this matches the way the BGEN file will be read in. 

The following shows how to run an analysis with a BGEN file using the dosage information. 


```julia
ordinalgwas(@formula(y ~ sex), datadir * "bgen_ex.csv", 
    datadir * "bgen_test"; geneticformat = "BGEN")
```




    StatsModels.TableRegressionModel{OrdinalMultinomialModel{Int64, Float64, LogitLink}, Matrix{Float64}}
    
    y ~ sex
    
    Coefficients:
    ──────────────────────────────────────────────────────
                   Estimate  Std.Error   t value  Pr(>|t|)
    ──────────────────────────────────────────────────────
    intercept1|2  -1.08828    0.128371  -8.47759    <1e-15
    intercept2|3   0.434834   0.118541   3.66821    0.0003
    intercept3|4   1.80196    0.145367  12.396      <1e-30
    sex            0.461387   0.162676   2.83623    0.0048
    ──────────────────────────────────────────────────────




```julia
run(`head ordinalgwas.pval.txt`);
```

    chr,pos,snpid,varid,maf,hwepval,infoscore,pval
    01,1001,RSID_101,SNPID_101,0.4169765,0.8627403,0.9846387832074154,0.12449779722495212
    01,2000,RSID_2,SNPID_2,0.19750981,0.1921813,9.0,0.0005572706473664106
    01,2001,RSID_102,SNPID_102,0.19766665,0.18440013,0.72730855295163,0.0004996980174197089
    01,3000,RSID_3,SNPID_3,0.48339605,0.9653543,0.955355323674357,0.5866179557677308
    01,3001,RSID_103,SNPID_103,0.4833961,0.9653536,0.9553553468765227,0.5866179209956521
    01,4000,RSID_4,SNPID_4,0.21670982,0.37192723,0.9917677420119885,0.3522822949083395
    01,4001,RSID_104,SNPID_104,0.2167098,0.37192774,0.99176792776617,0.35228235843382233
    01,5000,RSID_5,SNPID_5,0.38808233,0.5870128,0.9682577646241436,0.8349949433297409
    01,5001,RSID_105,SNPID_105,0.3880863,0.5867037,0.9682302221772928,0.8347201604324502


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




    StatsModels.TableRegressionModel{OrdinalMultinomialModel{Int64, Float64, ProbitLink}, Matrix{Float64}}
    
    trait ~ sex
    
    Coefficients:
    ──────────────────────────────────────────────────────
                   Estimate  Std.Error   t value  Pr(>|t|)
    ──────────────────────────────────────────────────────
    intercept1|2  -0.866156   0.210677  -4.11129    <1e-04
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
    1,758311,rs12562034,0.07763975155279501,0.4098763332666681,0.010076916742307648
    1,967643,rs2710875,0.32407407407407407,4.076249100705747e-7,2.62725649418621e-5
    1,1168108,rs11260566,0.19158878504672894,0.1285682279446898,1.0897484851087357e-5
    1,1375074,rs1312568,0.441358024691358,2.5376019650614977e-19,0.0051028839904383545
    1,1588771,rs35154105,0.0,1.0,1.0
    1,1789051,rs16824508,0.00462962962962965,0.9332783156468178,0.48653776297937934
    1,1990452,rs2678939,0.4537037037037037,5.07695957708431e-11,0.33231290090456195
    1,2194615,rs7553178,0.22685185185185186,0.17056143157457776,0.25915513977199217



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

In practice, we often perform GWAS on selected SNPs and/or selected samples. They can be specified by the `snpinds`, `covrowinds` and `geneticrowinds` keywords of `ordinalgwas` function. 

For example, to perform GWAS on SNPs with minor allele frequency (MAF) above 0.05


```julia
# create SNP mask
snpinds = maf(SnpArray(datadir * "hapmap3.bed")) .≥ 0.05
# GWAS on selected SNPs
@time ordinalgwas(@formula(trait ~ sex), datadir * "covariate.txt", datadir * "hapmap3", 
    snpinds=snpinds, nullfile="commonvariant.null.txt", pvalfile="commonvariant.pval.txt")
```

      1.578096 seconds (3.07 M allocations: 198.051 MiB, 8.05% gc time, 86.59% compilation time)





    StatsModels.TableRegressionModel{OrdinalMultinomialModel{Int64, Float64, LogitLink}, Matrix{Float64}}
    
    trait ~ sex
    
    Coefficients:
    ──────────────────────────────────────────────────────
                   Estimate  Std.Error   t value  Pr(>|t|)
    ──────────────────────────────────────────────────────
    intercept1|2  -1.48564    0.358891  -4.13952    <1e-04
    intercept2|3  -0.569479   0.341044  -1.66981    0.0959
    intercept3|4   0.429815   0.339642   1.26549    0.2066
    sex            0.424656   0.213914   1.98517    0.0480
    ──────────────────────────────────────────────────────




```julia
run(`head commonvariant.pval.txt`);
```

    chr,pos,snpid,maf,hwepval,pval
    1,758311,rs12562034,0.07763975155279501,0.4098763332666681,0.004565312839535881
    1,967643,rs2710875,0.32407407407407407,4.076249100705747e-7,3.108283828553597e-5
    1,1168108,rs11260566,0.19158878504672894,0.1285682279446898,1.2168672367654906e-5
    1,1375074,rs1312568,0.441358024691358,2.5376019650614977e-19,0.008206860046174836
    1,1990452,rs2678939,0.4537037037037037,5.07695957708431e-11,0.29972829571847304
    1,2194615,rs7553178,0.22685185185185186,0.17056143157457776,0.17133312458049288
    1,2396747,rs13376356,0.1448598130841121,0.9053079215078139,0.5320416198875104
    1,2823603,rs1563468,0.4830246913580247,4.23065537243926e-9,0.22519139178356798
    1,3025087,rs6690373,0.2538699690402477,9.238641887192776e-8,0.7018469417717436



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

`covrowinds` specify the samples in the covariate file and `geneticrowinds` for PLINK or VCF File. User should be particularly careful when using these two keyword. Selected rows in SnpArray should exactly match the samples in the null model. Otherwise the results are meaningless.

## Likelihood ratio test (LRT)

By default, `ordinalgwas` calculates p-value for each SNP using score test. Score test is fast because it doesn't require fitting alternative model for each SNP. User can request likelihood ratio test (LRT) using keyword `test=:lrt`. LRT is much slower but may be more powerful than score test.

Because LRT fits the alternative model for each SNP, we also output the standard error and estimated effect size of the SNP. 


```julia
@time ordinalgwas(@formula(trait ~ sex), datadir * "covariate.txt", datadir * "hapmap3", 
    test=:LRT, nullfile="lrt.null.txt", pvalfile="lrt.pval.txt")
```

     21.716944 seconds (5.51 M allocations: 1.911 GiB, 1.72% gc time, 1.84% compilation time)





    StatsModels.TableRegressionModel{OrdinalMultinomialModel{Int64, Float64, LogitLink}, Matrix{Float64}}
    
    trait ~ sex
    
    Coefficients:
    ──────────────────────────────────────────────────────
                   Estimate  Std.Error   t value  Pr(>|t|)
    ──────────────────────────────────────────────────────
    intercept1|2  -1.48564    0.358891  -4.13952    <1e-04
    intercept2|3  -0.569479   0.341044  -1.66981    0.0959
    intercept3|4   0.429815   0.339642   1.26549    0.2066
    sex            0.424656   0.213914   1.98517    0.0480
    ──────────────────────────────────────────────────────



Note the extra `effect` and `stder` columns in pvalfile, which is the effect size (regression coefficient) and standard errors for each SNP. `-1` is reported when there is no variation in the SNP (0 MAF).  


```julia
run(`head lrt.pval.txt`);
```

    chr,pos,snpid,maf,hwepval,effect,stder,pval
    1,554484,rs10458597,0.0,1.0,0.0,-1.0,1.0
    1,758311,rs12562034,0.07763975155279501,0.4098763332666681,-1.0057833719544336,0.34080917686608586,0.0019185836579804134
    1,967643,rs2710875,0.32407407407407407,4.076249100705747e-7,-0.6488560566291872,0.15659100275284915,1.805050556976133e-5
    1,1168108,rs11260566,0.19158878504672894,0.1285682279446898,-0.9157225669357891,0.2177988861470356,5.873384712686268e-6
    1,1375074,rs1312568,0.441358024691358,2.5376019650614977e-19,-0.3318136652577256,0.12697404813031218,0.008081022577833346
    1,1588771,rs35154105,0.0,1.0,0.0,-1.0,1.0
    1,1789051,rs16824508,0.00462962962962965,0.9332783156468178,-0.733802638869998,1.2495287571216083,0.5169027130144458
    1,1990452,rs2678939,0.4537037037037037,5.07695957708431e-11,-0.13586499231819757,0.13123472232496708,0.29946402200912603
    1,2194615,rs7553178,0.22685185185185186,0.17056143157457776,-0.251207564044013,0.17872139152979566,0.1615106909443865



```julia
# clean up
rm("lrt.pval.txt", force=true)
rm("lrt.null.txt", force=true)
```

In this example, GWAS by score test takes less than 0.2 second, while GWAS by LRT takes over 20 seconds. Over 100 fold difference in run time. 

## Score test for screening, LRT for power 

For large data sets, a practical solution is to perform the score test first across all SNPs, then re-do LRT for the most promising SNPs according to score test p-values.

**Step 1**: Perform score test GWAS, results in `score.pval.txt`.


```julia
@time ordinalgwas(@formula(trait ~ sex), datadir * "covariate.txt", datadir * "hapmap3", 
    test=:score, pvalfile="score.pval.txt");
```

      0.328529 seconds (500.80 k allocations: 49.171 MiB)



```julia
run(`head score.pval.txt`);
```

    chr,pos,snpid,maf,hwepval,pval
    1,554484,rs10458597,0.0,1.0,1.0
    1,758311,rs12562034,0.07763975155279501,0.4098763332666681,0.004565312839535881
    1,967643,rs2710875,0.32407407407407407,4.076249100705747e-7,3.108283828553597e-5
    1,1168108,rs11260566,0.19158878504672894,0.1285682279446898,1.2168672367654906e-5
    1,1375074,rs1312568,0.441358024691358,2.5376019650614977e-19,0.008206860046174836
    1,1588771,rs35154105,0.0,1.0,1.0
    1,1789051,rs16824508,0.00462962962962965,0.9332783156468178,0.5111981332534471
    1,1990452,rs2678939,0.4537037037037037,5.07695957708431e-11,0.29972829571847304
    1,2194615,rs7553178,0.22685185185185186,0.17056143157457776,0.17133312458049288


**Step 2**: Sort score test p-values and find top 10 SNPs.


```julia
scorepvals = CSV.read("score.pval.txt", DataFrame)[!, end] # p-values in last column
tophits = sortperm(scorepvals)[1:10] # indices of 10 SNPs with smallest p-values
scorepvals[tophits] # smallest 10 p-values
```




    10-element Vector{Float64}:
     1.3080149099173788e-6
     6.536722765046125e-6
     9.664742185663502e-6
     1.2168672367654906e-5
     1.8027460018322353e-5
     2.0989542284210004e-5
     2.6844521269932935e-5
     3.108283828553597e-5
     4.101091287507263e-5
     4.2966265138340985e-5



**Step 3**: Re-do LRT on top hits.


```julia
@time ordinalgwas(@formula(trait ~ sex), datadir * "covariate.txt", datadir * "hapmap3", 
    snpinds=tophits, test=:LRT, pvalfile="lrt.pval.txt");
```

      1.321790 seconds (1.91 M allocations: 114.008 MiB, 3.17% gc time, 95.63% compilation time)



```julia
run(`cat lrt.pval.txt`);
```

    chr,pos,snpid,maf,hwepval,effect,stder,pval
    1,967643,rs2710875,0.32407407407407407,4.076249100705747e-7,-0.6488560566291872,0.15659100275284915,1.805050556976133e-5
    1,1168108,rs11260566,0.19158878504672894,0.1285682279446898,-0.9157225669357891,0.2177988861470356,5.873384712686268e-6
    3,36821790,rs4678553,0.23456790123456794,0.1094668216324497,0.7424952268973514,0.17055414498937366,1.1303825016262592e-5
    4,11017683,rs16881446,0.27554179566563464,0.8942746118760274,-0.7870581482955514,0.18729099511224528,1.1105427468799613e-5
    5,3739190,rs12521166,0.0679012345679012,0.18613647746093887,1.1468852997925314,0.2809307861547482,4.7812882296576845e-5
    6,7574576,rs1885466,0.17746913580246915,0.7620687178987191,0.8750621092263018,0.19417734621377952,7.272346896740631e-6
    6,52474721,rs2073183,0.1826625386996904,0.5077765730476698,0.7790794914858653,0.19764574359737597,5.069394513904594e-5
    7,41152376,rs28880,0.3379629629629629,0.8052368892744096,-0.8146339024453504,0.17471459882610024,9.180126530294943e-7
    7,84223996,rs4128623,0.07870370370370372,0.0218347173467568,1.0022229316338562,0.25534365673199283,6.587895464657107e-5
    23,121048059,rs1937165,0.4380804953560371,3.959609737265113e-16,0.5392313636256605,0.1280959604207797,1.975464385552384e-5



```julia
# clean up
rm("ordinalgwas.null.txt", force=true)
rm("score.pval.txt", force=true)
rm("lrt.pval.txt", force=true)
```

## GxE or other interactions

### Testing jointly G + GxE 

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
    1,758311,rs12562034,0.07763975155279501,0.4098763332666681,0.017446010412245565
    1,967643,rs2710875,0.32407407407407407,4.076249100705747e-7,0.0001667073239488232
    1,1168108,rs11260566,0.19158878504672894,0.1285682279446898,4.7637624578926465e-5
    1,1375074,rs1312568,0.441358024691358,2.5376019650614977e-19,0.029138471242991387
    1,1588771,rs35154105,0.0,1.0,1.0
    1,1789051,rs16824508,0.00462962962962965,0.9332783156468178,0.29643631147339605
    1,1990452,rs2678939,0.4537037037037037,5.07695957708431e-11,0.3792458047934272
    1,2194615,rs7553178,0.22685185185185186,0.17056143157457776,0.32558226993239



```julia
# clean up
rm("ordinalgwas.null.txt")
rm("GxE.pval.txt")
```

### Testing only GxE interaction term

For some applications, the user may want to simply test the GxE interaction effect. This requires fitting the SNP in the null model and is much slower, but the command `ordinalgwas()` with keyword `analysistype = "gxe"` can be used test the interaction effect.
The environmental variable must be specified in the command using the keyword argument `e`, either as a symbol, such as `:age` or as a string `"age"`. 


```julia
ordinalgwas(@formula(trait ~ sex), datadir * "covariate.txt", datadir * "hapmap3"; analysistype = "gxe",
    e = :sex, pvalfile = "gxe_snp.pval.txt", snpinds=1:5, test=:score)
```




    StatsModels.TableRegressionModel{OrdinalMultinomialModel{Int64, Float64, LogitLink}, Matrix{Float64}}
    
    trait ~ sex
    
    Coefficients:
    ──────────────────────────────────────────────────────
                   Estimate  Std.Error   t value  Pr(>|t|)
    ──────────────────────────────────────────────────────
    intercept1|2  -1.48564    0.358891  -4.13952    <1e-04
    intercept2|3  -0.569479   0.341044  -1.66981    0.0959
    intercept3|4   0.429815   0.339642   1.26549    0.2066
    sex            0.424656   0.213914   1.98517    0.0480
    ──────────────────────────────────────────────────────




```julia
run(`head gxe_snp.pval.txt`);
```

    chr,pos,snpid,maf,hwepval,snpeffectnull,pval
    1,554484,rs10458597,0.0,1.0,0.0,1.0
    1,758311,rs12562034,0.07763975155279501,0.4098763332666681,-1.0057833719544336,0.6377422425980503
    1,967643,rs2710875,0.32407407407407407,4.076249100705747e-7,-0.6488560566291872,0.9667114197304937
    1,1168108,rs11260566,0.19158878504672894,0.1285682279446898,-0.9157225669357891,0.263526746940882
    1,1375074,rs1312568,0.441358024691358,2.5376019650614977e-19,-0.3318136652577256,0.7811133315583378



```julia
# clean up
rm("gxe_snp.pval.txt", force=true)
```

## SNP-set testing

In many applications, we want to test a SNP-set. The function with keyword `analysistype = "snpset"` can be used to do this. To specify the type of snpset test, use the `snpset` keyword argument. 

The snpset can be specified as either:
- a window (test every X snps) => `snpset = X`
- an annotated file.  This requires `snpset = filename` where filename is an input file, with no header and two columns separated by a space. The first column must contain the snpset ID and the second column must contain the snpid's (rsid) identical to the bimfile, or in the case of a BGEN format, the order the data will be read (see BGEN above for more details).
- a joint test on only a specific set of snps. `snpset = AbstractVector` where the vector specifies the snps you want to perform one joint snpset test for. The vector can either be a vector of integers where the elements are the indicies of the SNPs to test, a vector of booleans, where true represents that you want to select that SNP index, or a range indicating the indicies of the SNPs to test. 

In following example, we perform a SNP-set test on the 50th to 55th snps. 


```julia
ordinalgwas(@formula(trait ~ sex), datadir * "covariate.txt", datadir * "hapmap3"; 
    analysistype = "snpset", pvalfile = "snpset.pval.txt", snpset = 50:55)
```




    StatsModels.TableRegressionModel{OrdinalMultinomialModel{Int64, Float64, LogitLink}, Matrix{Float64}}
    
    trait ~ sex
    
    Coefficients:
    ──────────────────────────────────────────────────────
                   Estimate  Std.Error   t value  Pr(>|t|)
    ──────────────────────────────────────────────────────
    intercept1|2  -1.48564    0.358891  -4.13952    <1e-04
    intercept2|3  -0.569479   0.341044  -1.66981    0.0959
    intercept3|4   0.429815   0.339642   1.26549    0.2066
    sex            0.424656   0.213914   1.98517    0.0480
    ──────────────────────────────────────────────────────




```julia
run(`head snpset.pval.txt`);
```

    The joint pvalue of snps indexed at 50:55 is 0.3647126536663968



```julia
# clean up
rm("snpset.pval.txt", force=true)
rm("ordinalgwas.null.txt", force=true)
```

In the following example we run a SNP-set test on the annotated SNP-set file.


```julia
run(`head $datadir/hapmap_snpsetfile.txt`);
```

    gene1 rs10458597
    gene1 rs12562034
    gene1 rs2710875
    gene1 rs11260566
    gene1 rs1312568
    gene1 rs35154105
    gene1 rs16824508
    gene1 rs2678939
    gene1 rs7553178
    gene1 rs13376356



```julia
ordinalgwas(@formula(trait ~ sex), datadir * "covariate.txt", datadir * "hapmap3";
    analysistype = "snpset", pvalfile = "snpset.pval.txt", 
    snpset = datadir * "/hapmap_snpsetfile.txt")
```




    StatsModels.TableRegressionModel{OrdinalMultinomialModel{Int64, Float64, LogitLink}, Matrix{Float64}}
    
    trait ~ sex
    
    Coefficients:
    ──────────────────────────────────────────────────────
                   Estimate  Std.Error   t value  Pr(>|t|)
    ──────────────────────────────────────────────────────
    intercept1|2  -1.48564    0.358891  -4.13952    <1e-04
    intercept2|3  -0.569479   0.341044  -1.66981    0.0959
    intercept3|4   0.429815   0.339642   1.26549    0.2066
    sex            0.424656   0.213914   1.98517    0.0480
    ──────────────────────────────────────────────────────




```julia
run(`head snpset.pval.txt`);
```

    snpsetid,nsnps,pval
    gene1,93,1.7213394698790313e-5
    gene2,93,0.03692497684163047
    gene3,93,0.7478549371392569
    gene4,92,0.02765079822347653
    gene5,93,0.6119581594573154
    gene6,93,0.029184642230066303
    gene7,93,0.39160073488514746
    gene8,93,0.12539574109853385
    gene9,93,0.7085635621607591



```julia
# clean up
rm("snpset.pval.txt", force=true)
rm("ordinalgwas.null.txt", force=true)
```

In the following example we run a SNP-set test on every 15 SNPs.


```julia
ordinalgwas(@formula(trait ~ sex), datadir * "covariate.txt", datadir * "hapmap3";
    analysistype = "snpset", pvalfile = "snpset.pval.txt", snpset=15)
```




    StatsModels.TableRegressionModel{OrdinalMultinomialModel{Int64, Float64, LogitLink}, Matrix{Float64}}
    
    trait ~ sex
    
    Coefficients:
    ──────────────────────────────────────────────────────
                   Estimate  Std.Error   t value  Pr(>|t|)
    ──────────────────────────────────────────────────────
    intercept1|2  -1.48564    0.358891  -4.13952    <1e-04
    intercept2|3  -0.569479   0.341044  -1.66981    0.0959
    intercept3|4   0.429815   0.339642   1.26549    0.2066
    sex            0.424656   0.213914   1.98517    0.0480
    ──────────────────────────────────────────────────────




```julia
run(`head snpset.pval.txt`);
```

    startchr,startpos,startsnpid,endchr,endpos,endsnpid,pval
    1,554484,rs10458597,1,3431124,rs12093117,1.9394189465428366e-13
    1,3633945,rs10910017,1,6514524,rs932112,0.11077869162800988
    1,6715827,rs441515,1,9534606,rs4926480,0.2742450817958276
    1,9737551,rs12047054,1,12559747,rs4845907,0.49346113650467543
    1,12760427,rs848577,1,16021797,rs6679870,0.15447358658256338
    1,16228774,rs1972359,1,19100349,rs9426794,0.1544232917023417
    1,19301516,rs4912046,1,22122176,rs9426785,0.4749873349462694
    1,22323074,rs2235529,1,25166528,rs4648890,0.41246096621461137
    1,25368553,rs7527379,1,28208168,rs12140070,0.16458135278661626



```julia
# clean up
rm("snpset.pval.txt", force=true)
rm("ordinalgwas.null.txt", force=true)
```

## Plotting Results

To plot the GWAS results, we recommend using the [MendelPlots.jl package](https://openmendel.github.io/MendelPlots.jl/latest/).

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




    75-element Vector{String}:
     "/Users/xyz/.julia/dev/OrdinalGWAS/data/hapmap3.chr.1.bed"
     "/Users/xyz/.julia/dev/OrdinalGWAS/data/hapmap3.chr.1.bim"
     "/Users/xyz/.julia/dev/OrdinalGWAS/data/hapmap3.chr.1.fam"
     "/Users/xyz/.julia/dev/OrdinalGWAS/data/hapmap3.chr.10.bed"
     "/Users/xyz/.julia/dev/OrdinalGWAS/data/hapmap3.chr.10.bim"
     "/Users/xyz/.julia/dev/OrdinalGWAS/data/hapmap3.chr.10.fam"
     "/Users/xyz/.julia/dev/OrdinalGWAS/data/hapmap3.chr.11.bed"
     "/Users/xyz/.julia/dev/OrdinalGWAS/data/hapmap3.chr.11.bim"
     "/Users/xyz/.julia/dev/OrdinalGWAS/data/hapmap3.chr.11.fam"
     "/Users/xyz/.julia/dev/OrdinalGWAS/data/hapmap3.chr.12.bed"
     "/Users/xyz/.julia/dev/OrdinalGWAS/data/hapmap3.chr.12.bim"
     "/Users/xyz/.julia/dev/OrdinalGWAS/data/hapmap3.chr.12.fam"
     "/Users/xyz/.julia/dev/OrdinalGWAS/data/hapmap3.chr.13.bed"
     ⋮
     "/Users/xyz/.julia/dev/OrdinalGWAS/data/hapmap3.chr.6.bed"
     "/Users/xyz/.julia/dev/OrdinalGWAS/data/hapmap3.chr.6.bim"
     "/Users/xyz/.julia/dev/OrdinalGWAS/data/hapmap3.chr.6.fam"
     "/Users/xyz/.julia/dev/OrdinalGWAS/data/hapmap3.chr.7.bed"
     "/Users/xyz/.julia/dev/OrdinalGWAS/data/hapmap3.chr.7.bim"
     "/Users/xyz/.julia/dev/OrdinalGWAS/data/hapmap3.chr.7.fam"
     "/Users/xyz/.julia/dev/OrdinalGWAS/data/hapmap3.chr.8.bed"
     "/Users/xyz/.julia/dev/OrdinalGWAS/data/hapmap3.chr.8.bim"
     "/Users/xyz/.julia/dev/OrdinalGWAS/data/hapmap3.chr.8.fam"
     "/Users/xyz/.julia/dev/OrdinalGWAS/data/hapmap3.chr.9.bed"
     "/Users/xyz/.julia/dev/OrdinalGWAS/data/hapmap3.chr.9.bim"
     "/Users/xyz/.julia/dev/OrdinalGWAS/data/hapmap3.chr.9.fam"



Step 1: Fit the null model. Setting third argument `geneticfile` to `nothing` instructs `ordinalgwas` function to fit the null model only.


```julia
nm = ordinalgwas(@formula(trait ~ sex), datadir * "covariate.txt", nothing)
```




    StatsModels.TableRegressionModel{OrdinalMultinomialModel{Int64, Float64, LogitLink}, Matrix{Float64}}
    
    trait ~ sex
    
    Coefficients:
    ──────────────────────────────────────────────────────
                   Estimate  Std.Error   t value  Pr(>|t|)
    ──────────────────────────────────────────────────────
    intercept1|2  -1.48564    0.358891  -4.13952    <1e-04
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




    23-element Vector{String}:
     "/Users/xyz/.julia/dev/OrdinalGWAS/data/hapmap3.chr.1.pval.txt"
     "/Users/xyz/.julia/dev/OrdinalGWAS/data/hapmap3.chr.10.pval.txt"
     "/Users/xyz/.julia/dev/OrdinalGWAS/data/hapmap3.chr.11.pval.txt"
     "/Users/xyz/.julia/dev/OrdinalGWAS/data/hapmap3.chr.12.pval.txt"
     "/Users/xyz/.julia/dev/OrdinalGWAS/data/hapmap3.chr.13.pval.txt"
     "/Users/xyz/.julia/dev/OrdinalGWAS/data/hapmap3.chr.14.pval.txt"
     "/Users/xyz/.julia/dev/OrdinalGWAS/data/hapmap3.chr.15.pval.txt"
     "/Users/xyz/.julia/dev/OrdinalGWAS/data/hapmap3.chr.16.pval.txt"
     "/Users/xyz/.julia/dev/OrdinalGWAS/data/hapmap3.chr.17.pval.txt"
     "/Users/xyz/.julia/dev/OrdinalGWAS/data/hapmap3.chr.18.pval.txt"
     "/Users/xyz/.julia/dev/OrdinalGWAS/data/hapmap3.chr.19.pval.txt"
     "/Users/xyz/.julia/dev/OrdinalGWAS/data/hapmap3.chr.2.pval.txt"
     "/Users/xyz/.julia/dev/OrdinalGWAS/data/hapmap3.chr.20.pval.txt"
     "/Users/xyz/.julia/dev/OrdinalGWAS/data/hapmap3.chr.21.pval.txt"
     "/Users/xyz/.julia/dev/OrdinalGWAS/data/hapmap3.chr.22.pval.txt"
     "/Users/xyz/.julia/dev/OrdinalGWAS/data/hapmap3.chr.23.pval.txt"
     "/Users/xyz/.julia/dev/OrdinalGWAS/data/hapmap3.chr.3.pval.txt"
     "/Users/xyz/.julia/dev/OrdinalGWAS/data/hapmap3.chr.4.pval.txt"
     "/Users/xyz/.julia/dev/OrdinalGWAS/data/hapmap3.chr.5.pval.txt"
     "/Users/xyz/.julia/dev/OrdinalGWAS/data/hapmap3.chr.6.pval.txt"
     "/Users/xyz/.julia/dev/OrdinalGWAS/data/hapmap3.chr.7.pval.txt"
     "/Users/xyz/.julia/dev/OrdinalGWAS/data/hapmap3.chr.8.pval.txt"
     "/Users/xyz/.julia/dev/OrdinalGWAS/data/hapmap3.chr.9.pval.txt"



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
run(`cat $datadir/../docs/cluster_preparedata.jl`);
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
run(`cat $datadir/../docs/cluster_run.jl`);
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


## Proportional Odds Assumption

The ordered multinomial model used in this package is a cumulative logit model where one effect size is estimated for each covariate. Using the logit link yields the commonly used proportional odds model which assumes the proportional odds assumption holds for your data/model. The proportional odds assumption assumes that the effect of a covariate $\boldsymbol{\beta}$ is constant across different groupings of the ordinal outcome i.e.

The SNP's effect (on the odds ratio) is the same for the odds ratio of mild vs {medium, severe} as it is for the odds ratio of {mild, medium} vs severe. 

The consequences/neccesity of this assumption has been thoroughly discussed, as an example see [this blog post](https://www.fharrell.com/post/po/). In simualtions, we did not find violation of this to increase type I error. As stated in the post, when the assumption is violated (meaning there are different effect sizes for different groupings), the estimated effect size is like a weighted average of the different effect sizes. Formal tests for proportional odds are likely to reject the assumption for large n (see p. 306, Categorical Data Analysis 3rd edition by Agresti). Violation of this assumption means the effect size $\boldsymbol{\beta}$ differs across outcome strata, but should not imply a significant finding is not significant.
