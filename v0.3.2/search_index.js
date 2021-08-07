var documenterSearchIndex = {"docs": [

{
    "location": "#",
    "page": "OrdinalGWAS.jl",
    "title": "OrdinalGWAS.jl",
    "category": "page",
    "text": ""
},

{
    "location": "#OrdinalGWAS.jl-1",
    "page": "OrdinalGWAS.jl",
    "title": "OrdinalGWAS.jl",
    "category": "section",
    "text": "OrdinalGWAS.jl is a Julia package for performing genome-wide association studies (GWAS) for ordered categorical phenotypes using proportional odds model or ordred Probit model. It is useful when the phenotype takes ordered discrete values, e.g., disease status (undiagnosed, pre-disease, mild, moderate, severe)."
},

{
    "location": "#Installation-1",
    "page": "OrdinalGWAS.jl",
    "title": "Installation",
    "category": "section",
    "text": "This package requires Julia v0.7 or later and two other unregistered packages SnpArrays.jl and OrdinalMultinomialModels.jl. The package has not yet been registered and must be installed using the repository location. Start julia and use the ] key to switch to the package manager REPL(v1.0) pkg> add https://github.com/OpenMendel/SnpArrays.jl.git\n(v1.0) pkg> add https://github.com/OpenMendel/OrdinalMultinomialModels.jl.git\n(v1.0) pkg> add https://github.com/OpenMendel/OrdinalGWAS.jl.git# machine information for this tutorial\nversioninfo()Julia Version 1.0.3\nCommit 099e826241 (2018-12-18 01:34 UTC)\nPlatform Info:\n  OS: macOS (x86_64-apple-darwin14.5.0)\n  CPU: Intel(R) Core(TM) i7-6920HQ CPU @ 2.90GHz\n  WORD_SIZE: 64\n  LIBM: libopenlibm\n  LLVM: libLLVM-6.0.0 (ORCJIT, skylake)\nEnvironment:\n  JULIA_EDITOR = code# for use in this tutorial\nusing BenchmarkTools, CSV, Glob, OrdinalGWAS, SnpArrays"
},

{
    "location": "#Example-data-set-1",
    "page": "OrdinalGWAS.jl",
    "title": "Example data set",
    "category": "section",
    "text": "data folder of the package contains an example data set. In general, user can locate this folder by commandusing OrdinalGWAS\nconst datadir = normpath(joinpath(dirname(pathof(OrdinalGWAS)), \"../data/\"))\"/Users/huazhou/.julia/dev/OrdinalGWAS.jl/data/\"# content of the data folder\nreaddir(glob\"*.*\", datadir)6-element Array{String,1}:\n \"/Users/huazhou/.julia/dev/OrdinalGWAS.jl/data/covariate.txt\"\n \"/Users/huazhou/.julia/dev/OrdinalGWAS.jl/data/hapmap3.bed\"  \n \"/Users/huazhou/.julia/dev/OrdinalGWAS.jl/data/hapmap3.bim\"  \n \"/Users/huazhou/.julia/dev/OrdinalGWAS.jl/data/hapmap3.fam\"  \n \"/Users/huazhou/.julia/dev/OrdinalGWAS.jl/data/hapmap3.map\"  \n \"/Users/huazhou/.julia/dev/OrdinalGWAS.jl/data/simtrait.jl\""
},

{
    "location": "#OrdinalGWAS.ordinalgwas",
    "page": "OrdinalGWAS.jl",
    "title": "OrdinalGWAS.ordinalgwas",
    "category": "function",
    "text": "ordinalgwas(nullformula, covfile, plinkfile)\nordinalgwas(nullformula, df, plinkfile)\nordinalgwas(fittednullmodel, plinkfile)\nordinalgwas(fittednullmodel, bedfile, bimfile, bedn)\n\nPositional arguments\n\nnullformula::Formula: formula for the null model.\ncovfile::AbstractString: covariate file (csv) with one header line. One column    should be the ordered categorical phenotype coded as integers starting from 1.\ndf::DataFrame: DataFrame containing response and regressors for null model.\nplinkfile::AbstractString: Plink file name without the bed, fam, or bim    extensions. If plinkfile==nothing, only null model is fitted. If plinkfile    is provided, bed, bim, and fam file with same plinkfile prefix need to exist.    Compressed file formats such as gz and bz2 are allowed. Check all allowed formats    by SnpArrays.ALLOWED_FORMAT.  \nfittednullmodel::StatsModels.DataFrameRegressionModel: the fitted null model    output from ordinalgwas(nullformula, covfile) or ordinalgwas(nullformula, df).\nbedfile::Union{AbstractString,IOStream}: path to Plink bed file.\nbimfile::Union{AbstractString,IOStream}: path to Plink bim file.\nbedn::Integer: number of samples in bed file.\n\nKeyword arguments\n\nnullfile::Union{AbstractString, IOStream}: output file for the fitted null model; default is    ordinalgwas.null.txt. \npvalfile::Union{AbstractString, IOStream}: output file for the gwas p-values; default is    ordinalgwas.pval.txt. \ncovtype::Vector{DataType}: type information for covfile. This is useful   when CSV.read(covarfile) has parsing errors.  \ncovrowinds::Union{Nothing,AbstractVector{<:Integer}}: sample indices for covariate file.  \ntestformula::Formula: formula for test unit. Default is @formula(trait ~ 0 + snp).\ntest::Symbol: :score (default) or :lrt.  \nlink::GLM.Link: LogitLink() (default), ProbitLink(), CauchitLink(),   or CloglogLink().\nsnpmodel: ADDITIVE_MODEL (default), DOMINANT_MODEL, or RECESSIVE_MODEL.\nsnpinds::Union{Nothing,AbstractVector{<:Integer}}: SNP indices for bed file.\nbedrowinds::Union{Nothing,AbstractVector{<:Integer}}: sample indices for bed file.\nsolver: an optimization solver supported by MathProgBase. Default is    NLoptSolver(algorithm=:LD_SLSQP, maxeval=4000). Another common choice is    IpoptSolver(print_level=0).\nverbose::Bool: default is false.\n\n\n\n\n\n"
},

{
    "location": "#Basic-usage-1",
    "page": "OrdinalGWAS.jl",
    "title": "Basic usage",
    "category": "section",
    "text": "The following command performs GWAS using the proportional odds model. The output is the fitted null model.ordinalgwas(@formula(trait ~ sex), datadir * \"covariate.txt\", datadir * \"hapmap3\")StatsModels.DataFrameRegressionModel{OrdinalMultinomialModel{Int64,Float64,LogitLink},Array{Float64,2}}\n\nFormula: trait ~ +sex\n\nCoefficients:\n      Estimate Std.Error  t value Pr(>|t|)\nθ1    -1.48564  0.358891 -4.13952    <1e-4\nθ2   -0.569479  0.341044 -1.66981   0.0959\nθ3    0.429815  0.339642  1.26549   0.2066\nβ1    0.424656  0.213914  1.98517   0.0480For documentation of the ordinalgwas function, type ?ordinalgwas in Julia REPL.ordinalgwas"
},

{
    "location": "#Formula-for-null-model-1",
    "page": "OrdinalGWAS.jl",
    "title": "Formula for null model",
    "category": "section",
    "text": "The first argument specifies the null model without SNP effects, e.g., @formula(trait ~ sex)."
},

{
    "location": "#Input-files-1",
    "page": "OrdinalGWAS.jl",
    "title": "Input files",
    "category": "section",
    "text": "ordinalgwas expects two input files: one for responses plus covariates (second argument), the other the Plink files for genotypes (third argument)."
},

{
    "location": "#Covariate-and-trait-file-1",
    "page": "OrdinalGWAS.jl",
    "title": "Covariate and trait file",
    "category": "section",
    "text": "Covariates and phenotype are provided in a csv file, e.g., covariate.txt, which has one header line for variable names. In this example, variable trait is the ordered categorical phenotypes coded as integers 1 to 4. We want to include variable sex as the covariate in GWAS.run(`head $(datadir)covariate.txt`);famid,perid,faid,moid,sex,trait\n2431,NA19916,0,0,1,4\n2424,NA19835,0,0,2,4\n2469,NA20282,0,0,2,4\n2368,NA19703,0,0,1,3\n2425,NA19901,0,0,2,3\n2427,NA19908,0,0,1,4\n2430,NA19914,0,0,2,4\n2470,NA20287,0,0,2,1\n2436,NA19713,0,0,2,3"
},

{
    "location": "#Plink-file-1",
    "page": "OrdinalGWAS.jl",
    "title": "Plink file",
    "category": "section",
    "text": "Genotype data is available as binary Plink files.readdir(glob\"hapmap3.*\", datadir)4-element Array{String,1}:\n \"/Users/huazhou/.julia/dev/OrdinalGWAS.jl/data/hapmap3.bed\"\n \"/Users/huazhou/.julia/dev/OrdinalGWAS.jl/data/hapmap3.bim\"\n \"/Users/huazhou/.julia/dev/OrdinalGWAS.jl/data/hapmap3.fam\"\n \"/Users/huazhou/.julia/dev/OrdinalGWAS.jl/data/hapmap3.map\"In this example, there are 324 samples at 13,928 SNPs.size(SnpArray(datadir * \"hapmap3.bed\"))(324, 13928)Compressed Plink files are supported. For example, if Plink files are hapmap3.bed.gz, hapmap3.bim.gz and hapmap3.fam.gz, the same commandordinalgwas(@formula(trait ~ sex), datadir * \"covariate.txt\", datadir * \"hapmap3\")still works. Check all supported compression format bySnpArrays.ALLOWED_FORMAT6-element Array{String,1}:\n \"gz\"  \n \"zlib\"\n \"zz\"  \n \"xz\"  \n \"zst\" \n \"bz2\""
},

{
    "location": "#Output-files-1",
    "page": "OrdinalGWAS.jl",
    "title": "Output files",
    "category": "section",
    "text": "ordinalgwas outputs two files: ordinalgwas.null.txt and ordinalgwas.pval.txt. ordinalgwas.null.txt lists the estimated null model (without SNPs). run(`cat ordinalgwas.null.txt`);StatsModels.DataFrameRegressionModel{OrdinalMultinomialModel{Int64,Float64,LogitLink},Array{Float64,2}}\n\nFormula: trait ~ +sex\n\nCoefficients:\n      Estimate Std.Error  t value Pr(>|t|)\nθ1    -1.48564  0.358891 -4.13952    <1e-4\nθ2   -0.569479  0.341044 -1.66981   0.0959\nθ3    0.429815  0.339642  1.26549   0.2066\nβ1    0.424656  0.213914  1.98517   0.0480ordinalgwas.pval.txt tallies the SNPs and their pvalues. run(`head ordinalgwas.pval.txt`);chr,pos,snpid,maf,pval\n1,554484,rs10458597,0.0,1.0\n1,758311,rs12562034,0.07763975155279501,0.004565312839540994\n1,967643,rs2710875,0.32407407407407407,3.1082838285548695e-5\n1,1168108,rs11260566,0.19158878504672894,1.2168672367668912e-5\n1,1375074,rs1312568,0.441358024691358,0.008206860046175221\n1,1588771,rs35154105,0.0,1.0\n1,1789051,rs16824508,0.00462962962962965,0.5111981332544\n1,1990452,rs2678939,0.4537037037037037,0.29972829571847825\n1,2194615,rs7553178,0.22685185185185186,0.1713331245805063Output file names can be changed by the nullfile and pvalfile keywords respectively. For example, ordinalgwas(@formula(trait ~ sex), datadir * \"covariate.txt\", datadir * \"hapmap3\", pvalfile=\"ordinalgwas.pval.txt.gz\")will output the p-value file in compressed gz format."
},

{
    "location": "#Subsamples-1",
    "page": "OrdinalGWAS.jl",
    "title": "Subsamples",
    "category": "section",
    "text": "Use the keyword covrowinds to specify selected samples in the covarite file. Use the keyword bedrowinds to specify selected samples in the Plink bed file. For example, to use the first 300 samples in both covariate and bed file:ordinalgwas(@formula(trait ~ sex), covfile, plkfile, covrowinds=1:300, bedrowinds=1:300)note: Note\nUsers should always make sure that the selected samples in covariate file match exactly those in bed file. "
},

{
    "location": "#Input-non-genetic-data-as-DataFrame-1",
    "page": "OrdinalGWAS.jl",
    "title": "Input non-genetic data as DataFrame",
    "category": "section",
    "text": "Internally ordinalgwas parses the covariate file as a DataFrame by CSV.read(covfile). For covariate file of other formats, users can parse it as a DataFrame and then input the DataFrame to ordinalgwas directly.ordinalgwas(@formula(trait ~ sex), df, plinkfile)note: Note\nUsers should always make sure that individuals in covariate file or DataFrame match those in Plink fam file. For example, following code checks that the first 2 columns of the covariate.txt file match the first 2 columns of the hapmap3.fam file exactly.covdf = CSV.read(datadir * \"covariate.txt\")\nplkfam = CSV.read(datadir * \"hapmap3.fam\", header=0, delim=\' \')\nall(covdf[1] .== plkfam[1]) && all(covdf[2] .== plkfam[2])true"
},

{
    "location": "#Timing-1",
    "page": "OrdinalGWAS.jl",
    "title": "Timing",
    "category": "section",
    "text": "For this moderate-sized data set, ordinalgwas takes less than 0.2 second.@btime(ordinalgwas(@formula(trait ~ sex), datadir * \"covariate.txt\", datadir * \"hapmap3\"));  162.665 ms (710982 allocations: 33.35 MiB)# clean up\nrm(\"ordinalgwas.null.txt\")\nrm(\"ordinalgwas.pval.txt\")"
},

{
    "location": "#Link-functions-1",
    "page": "OrdinalGWAS.jl",
    "title": "Link functions",
    "category": "section",
    "text": "The link keyword argument of ordinalgwas can take value:  LogitLink(), proportional odds model (default),  \nProbitLink(), ordred Probit model,  \nCloglogLink(), proportional hazards model, or \nCauchyLink().For example, to perform GWAS using the ordred Probit modelordinalgwas(@formula(trait ~ sex), datadir * \"covariate.txt\", datadir * \"hapmap3\", \n    link=ProbitLink(), nullfile=\"opm.null.txt\", pvalfile=\"opm.pval.txt\")StatsModels.DataFrameRegressionModel{OrdinalMultinomialModel{Int64,Float64,ProbitLink},Array{Float64,2}}\n\nFormula: trait ~ +sex\n\nCoefficients:\n      Estimate Std.Error  t value Pr(>|t|)\nθ1   -0.866156  0.210677 -4.11129    <1e-4\nθ2   -0.359878  0.205817 -1.74854   0.0813\nθ3    0.247054  0.205382   1.2029   0.2299\nβ1    0.251058  0.128225  1.95795   0.0511The estimates in null model and p-values are slightly different from those in proportional odds moodel.run(`head opm.pval.txt`);chr,pos,snpid,maf,pval\n1,554484,rs10458597,0.0,1.0\n1,758311,rs12562034,0.07763975155279501,0.010076916742300138\n1,967643,rs2710875,0.32407407407407407,2.6272564941853933e-5\n1,1168108,rs11260566,0.19158878504672894,1.0897484851078458e-5\n1,1375074,rs1312568,0.441358024691358,0.005102883990438149\n1,1588771,rs35154105,0.0,1.0\n1,1789051,rs16824508,0.00462962962962965,0.48653776297859236\n1,1990452,rs2678939,0.4537037037037037,0.33231290090455434\n1,2194615,rs7553178,0.22685185185185186,0.25915513977197435rm(\"opm.null.txt\")\nrm(\"opm.pval.txt\")"
},

{
    "location": "#SNP-models-1",
    "page": "OrdinalGWAS.jl",
    "title": "SNP models",
    "category": "section",
    "text": "Genotypes are translated into numeric values according to different genetic model, which is specified by the snpmodel keyword. Default is ADDITIVE_MODEL.Genotype SnpArray ADDITIVE_MODEL DOMINANT_MODEL RECESSIVE_MODEL\nA1,A1 0x00 0 0 0\nmissing 0x01 NaN NaN NaN\nA1,A2 0x02 1 1 0\nA2,A2 0x03 2 1 1note: Note\nordinalgwas imputes missing genotypes according to minor allele frequencies. Users are advised to impute genotypes using more sophiscated methods before GWAS."
},

{
    "location": "#SNP-and/or-sample-masks-1",
    "page": "OrdinalGWAS.jl",
    "title": "SNP and/or sample masks",
    "category": "section",
    "text": "In practice, we often perform GWAS on selected SNPs and/or selected samples. They can be specified by the snpinds, covrowinds and bedrowinds keywords of ordinalgwas function. For example, to perform GWAS on SNPs with minor allele frequency (MAF) above 0.05# create SNP mask\nsnpinds = maf(SnpArray(\"../data/hapmap3.bed\")) .≥ 0.05\n# GWAS on selected SNPs\n@time ordinalgwas(@formula(trait ~ sex), datadir * \"covariate.txt\", datadir * \"hapmap3\", \n    snpinds=snpinds, nullfile=\"commonvariant.null.txt\", pvalfile=\"commonvariant.pval.txt\")  0.296101 seconds (881.81 k allocations: 42.526 MiB, 4.20% gc time)\n\n\n\n\n\nStatsModels.DataFrameRegressionModel{OrdinalMultinomialModel{Int64,Float64,LogitLink},Array{Float64,2}}\n\nFormula: trait ~ +sex\n\nCoefficients:\n      Estimate Std.Error  t value Pr(>|t|)\nθ1    -1.48564  0.358891 -4.13952    <1e-4\nθ2   -0.569479  0.341044 -1.66981   0.0959\nθ3    0.429815  0.339642  1.26549   0.2066\nβ1    0.424656  0.213914  1.98517   0.0480run(`head commonvariant.pval.txt`);chr,pos,snpid,maf,pval\n1,758311,rs12562034,0.07763975155279501,0.004565312839540994\n1,967643,rs2710875,0.32407407407407407,3.1082838285548695e-5\n1,1168108,rs11260566,0.19158878504672894,1.2168672367668912e-5\n1,1375074,rs1312568,0.441358024691358,0.008206860046175221\n1,1990452,rs2678939,0.4537037037037037,0.29972829571847825\n1,2194615,rs7553178,0.22685185185185186,0.1713331245805063\n1,2396747,rs13376356,0.1448598130841121,0.5320416198875456\n1,2823603,rs1563468,0.4830246913580247,0.225191391783573\n1,3025087,rs6690373,0.2538699690402477,0.7018469417717486# extra header line in commonvariant.pval.txt\ncountlines(\"commonvariant.pval.txt\"), count(snpinds)(12086, 12085)# clean up\nrm(\"commonvariant.null.txt\")\nrm(\"commonvariant.pval.txt\")covrowinds specify the samples in the covariate file and bedrowinds for SnpArray. User should be particularly careful when these two keyword. Selected rows in SnpArray should exactly match the samples in the null model. Otherwise the results are meaningless."
},

{
    "location": "#Likelihood-ratio-test-(LRT)-1",
    "page": "OrdinalGWAS.jl",
    "title": "Likelihood ratio test (LRT)",
    "category": "section",
    "text": "By default, ordinalgwas calculates p-value for each SNP using score test. Score test is fast because it doesn\'t require fitting alternative model for each SNP. User can request likelihood ratio test (LRT) using keyword test=:lrt. LRT is much slower but may be more powerful than score test.@time ordinalgwas(@formula(trait ~ sex), datadir * \"covariate.txt\", datadir * \"hapmap3\", \n    test=:LRT, nullfile=\"lrt.null.txt\", pvalfile=\"lrt.pval.txt\") 19.769862 seconds (8.18 M allocations: 2.044 GiB, 2.11% gc time)\n\n\n\n\n\nStatsModels.DataFrameRegressionModel{OrdinalMultinomialModel{Int64,Float64,LogitLink},Array{Float64,2}}\n\nFormula: trait ~ +sex\n\nCoefficients:\n      Estimate Std.Error  t value Pr(>|t|)\nθ1    -1.48564  0.358891 -4.13952    <1e-4\nθ2   -0.569479  0.341044 -1.66981   0.0959\nθ3    0.429815  0.339642  1.26549   0.2066\nβ1    0.424656  0.213914  1.98517   0.0480Note the extra effect column in pvalfile, which is the effect size (regression coefficient) for each SNP. run(`head lrt.pval.txt`);chr,pos,snpid,maf,effect,pval\n1,554484,rs10458597,0.0,0.0,1.0\n1,758311,rs12562034,0.07763975155279501,-1.0057833719544331,0.0019185836579804134\n1,967643,rs2710875,0.32407407407407407,-0.6488560566295055,1.805050556976241e-5\n1,1168108,rs11260566,0.19158878504672894,-0.9157225669357879,5.873384712685568e-6\n1,1375074,rs1312568,0.441358024691358,-0.33181366525772593,0.008081022577832324\n1,1588771,rs35154105,0.0,0.0,1.0\n1,1789051,rs16824508,0.00462962962962965,-0.7338026388701573,0.5169027130129711\n1,1990452,rs2678939,0.4537037037037037,-0.13586499231819726,0.29946402200912603\n1,2194615,rs7553178,0.22685185185185186,-0.2512075640440123,0.16151069094439868# clean up\nrm(\"lrt.pval.txt\")\nrm(\"lrt.null.txt\")In this example, GWAS by score test takes less than 0.2 second, while GWAS by LRT takes about 20 seconds. About 100 fold difference in run time. "
},

{
    "location": "#Score-test-for-screening,-LRT-for-power-1",
    "page": "OrdinalGWAS.jl",
    "title": "Score test for screening, LRT for power",
    "category": "section",
    "text": "For large data sets, a practical solution is to perform score test first, then re-do LRT for the most promising SNPs according to score test p-values.Step 1: Perform score test GWAS, results in score.pval.txt.@time ordinalgwas(@formula(trait ~ sex), datadir * \"covariate.txt\", datadir * \"hapmap3\", \n    test=:score, pvalfile=\"score.pval.txt\");  0.246964 seconds (758.61 k allocations: 35.808 MiB, 10.78% gc time)run(`head score.pval.txt`);chr,pos,snpid,maf,pval\n1,554484,rs10458597,0.0,1.0\n1,758311,rs12562034,0.07763975155279501,0.004565312839540994\n1,967643,rs2710875,0.32407407407407407,3.1082838285548695e-5\n1,1168108,rs11260566,0.19158878504672894,1.2168672367668912e-5\n1,1375074,rs1312568,0.441358024691358,0.008206860046175221\n1,1588771,rs35154105,0.0,1.0\n1,1789051,rs16824508,0.00462962962962965,0.5111981332544\n1,1990452,rs2678939,0.4537037037037037,0.29972829571847825\n1,2194615,rs7553178,0.22685185185185186,0.1713331245805063Step 2: Sort score test p-values and find top 10 SNPs.scorepvals = CSV.read(\"score.pval.txt\")[5] # p-values in 5th column\ntophits = sortperm(scorepvals)[1:10] # indices of 10 SNPs with smallest p-values\nscorepvals[tophits] # smallest 10 p-values10-element Array{Union{Missing, Float64},1}:\n 1.3080149099181335e-6\n 6.536722765052079e-6 \n 9.664742185669054e-6 \n 1.2168672367668912e-5\n 1.802746001833127e-5 \n 2.0989542284213636e-5\n 2.6844521269963608e-5\n 3.1082838285548695e-5\n 4.1010912875160476e-5\n 4.2966265138454806e-5Step 3: Re-do LRT on top hits.@time ordinalgwas(@formula(trait ~ sex), datadir * \"covariate.txt\", datadir * \"hapmap3\", \n    snpinds=tophits, test=:LRT, pvalfile=\"lrt.pval.txt\");  0.215692 seconds (358.45 k allocations: 20.114 MiB, 3.48% gc time)run(`cat lrt.pval.txt`);chr,pos,snpid,maf,effect,pval\n1,967643,rs2710875,0.32407407407407407,-0.6488560566295055,1.805050556976241e-5\n1,1168108,rs11260566,0.19158878504672894,-0.9157225669357879,5.873384712685568e-6\n3,36821790,rs4678553,0.23456790123456794,0.7424952268973518,1.1303825016262592e-5\n4,11017683,rs16881446,0.27554179566563464,-0.7870581482955515,1.1105427468799613e-5\n5,3739190,rs12521166,0.0679012345679012,1.1468852997925316,4.781288229657399e-5\n6,7574576,rs1885466,0.17746913580246915,0.8750621092263019,7.272346896740631e-6\n6,52474721,rs2073183,0.1826625386996904,0.7790794914858663,5.069394513906121e-5\n7,41152376,rs28880,0.3379629629629629,-0.814633902445351,9.180126530294943e-7\n7,84223996,rs4128623,0.07870370370370372,1.0022229316338573,6.587895464657512e-5\n23,121048059,rs1937165,0.4380804953560371,0.5392313636256612,1.9754643855522616e-5# clean up\nrm(\"ordinalgwas.null.txt\")\nrm(\"score.pval.txt\")\nrm(\"lrt.pval.txt\")"
},

{
    "location": "#GxE-or-other-interactions-1",
    "page": "OrdinalGWAS.jl",
    "title": "GxE or other interactions",
    "category": "section",
    "text": "In many applications, we want to test SNP effect and/or its interaction with other terms. testformula keyword specifies the test unit besides the covariates in nullformula. In following example, keyword testformula=@formula(trait ~ snp + snp & sex) instructs ordinalgwas to test joint effect of snp and snp & sex interaction.ordinalgwas(@formula(trait ~ sex), datadir * \"covariate.txt\", datadir * \"hapmap3\", \n    pvalfile=\"GxE.pval.txt\", testformula=@formula(trait ~ snp + snp & sex));run(`head GxE.pval.txt`);chr,pos,snpid,maf,pval\n1,554484,rs10458597,0.0,1.0\n1,758311,rs12562034,0.07763975155279501,0.017446010412254197\n1,967643,rs2710875,0.32407407407407407,0.0001667073239489097\n1,1168108,rs11260566,0.19158878504672894,4.763762457893366e-5\n1,1375074,rs1312568,0.441358024691358,0.029138471242993652\n1,1588771,rs35154105,0.0,1.0\n1,1789051,rs16824508,0.00462962962962965,0.2964363114944328\n1,1990452,rs2678939,0.4537037037037037,0.37924580479348785\n1,2194615,rs7553178,0.22685185185185186,0.325582269932396# clean up\nrm(\"ordinalgwas.null.txt\")\nrm(\"GxE.pval.txt\")"
},

{
    "location": "#Plotting-Results-1",
    "page": "OrdinalGWAS.jl",
    "title": "Plotting Results",
    "category": "section",
    "text": "To plot the GWAS results, use the MendelPlots.jl package."
},

{
    "location": "#Docker-1",
    "page": "OrdinalGWAS.jl",
    "title": "Docker",
    "category": "section",
    "text": "For ease of using OrdinalGWAS, we provide a Dockerfile so users don\'t need to install Julia and required packages. Only Docker app needs to be installed in order to run analysis. Following is tested on Docker 2.0.0.0-mac78.Step 1: Create a Dockerfile with content here, or, if the bash command wget is available, obtain Dockerfile by# on command line\nwget https://raw.githubusercontent.com/OpenMendel/OrdinalGWAS.jl/master/docker/DockerfileStep 2: Build a docker image called ordinalgwas-app, assuming that the Dockerfile is located in the ../docker folder. Building the image for the first time can take up to 10 minutes; but it only needs to be done once.# on command line\ndocker build -t ordinalgwas-app ../docker/Step 3: Suppose data files are located at /path/to/data folder, run analysis by# on command line\ndocker run -v /path/to/data:/data -t ordinalgwas-app julia -e \'using OrdinalGWAS; ordinalgwas(@formula(trait ~ sex), \"/data/covariate.txt\", \"/data/hapmap3\", nullfile=\"/data/ordinalgwas.null.txt\", pvalfile=\"/data/ordinalgwas\");\'Here  -t ordinalgwas-app creates a container using the ordinalgwas-app image build in step 2.  \n-v /path/to/data:/data maps the /path/to/data folder on host machine to the /data folder within the container. \njulia -e \'using OrdinalGWAS; ordinalgwas(@formula(trait ~ sex), \"/data/covariate.txt\", \"/data/hapmap3\", nullfile=\"/data/ordinalgwas.null.txt\", pvalfile=\"/data/ordinalgwas\"); calls Julia and runs ordinalgwas function. The output files are written in /path/to/data directory."
},

{
    "location": "#Multiple-Plink-file-sets-1",
    "page": "OrdinalGWAS.jl",
    "title": "Multiple Plink file sets",
    "category": "section",
    "text": "In large scale studies, genotypes data are split into multiple Plink files, e.g., by chromosome. Then GWAS analysis can be done in parallel. This can be achieved by two steps.Let\'s first create demo data by splitting hapmap3 according to chromosome:# split example hapmap3 data according to chromosome\nSnpArrays.split_plink(datadir * \"hapmap3\", :chromosome; prefix=datadir * \"hapmap3.chr.\")\nreaddir(glob\"hapmap3.chr.*\", datadir)75-element Array{String,1}:\n \"/Users/huazhou/.julia/dev/OrdinalGWAS.jl/data/hapmap3.chr.1.bed\" \n \"/Users/huazhou/.julia/dev/OrdinalGWAS.jl/data/hapmap3.chr.1.bim\" \n \"/Users/huazhou/.julia/dev/OrdinalGWAS.jl/data/hapmap3.chr.1.fam\" \n \"/Users/huazhou/.julia/dev/OrdinalGWAS.jl/data/hapmap3.chr.10.bed\"\n \"/Users/huazhou/.julia/dev/OrdinalGWAS.jl/data/hapmap3.chr.10.bim\"\n \"/Users/huazhou/.julia/dev/OrdinalGWAS.jl/data/hapmap3.chr.10.fam\"\n \"/Users/huazhou/.julia/dev/OrdinalGWAS.jl/data/hapmap3.chr.11.bed\"\n \"/Users/huazhou/.julia/dev/OrdinalGWAS.jl/data/hapmap3.chr.11.bim\"\n \"/Users/huazhou/.julia/dev/OrdinalGWAS.jl/data/hapmap3.chr.11.fam\"\n \"/Users/huazhou/.julia/dev/OrdinalGWAS.jl/data/hapmap3.chr.12.bed\"\n \"/Users/huazhou/.julia/dev/OrdinalGWAS.jl/data/hapmap3.chr.12.bim\"\n \"/Users/huazhou/.julia/dev/OrdinalGWAS.jl/data/hapmap3.chr.12.fam\"\n \"/Users/huazhou/.julia/dev/OrdinalGWAS.jl/data/hapmap3.chr.13.bed\"\n ⋮                                                                 \n \"/Users/huazhou/.julia/dev/OrdinalGWAS.jl/data/hapmap3.chr.6.bed\" \n \"/Users/huazhou/.julia/dev/OrdinalGWAS.jl/data/hapmap3.chr.6.bim\" \n \"/Users/huazhou/.julia/dev/OrdinalGWAS.jl/data/hapmap3.chr.6.fam\" \n \"/Users/huazhou/.julia/dev/OrdinalGWAS.jl/data/hapmap3.chr.7.bed\" \n \"/Users/huazhou/.julia/dev/OrdinalGWAS.jl/data/hapmap3.chr.7.bim\" \n \"/Users/huazhou/.julia/dev/OrdinalGWAS.jl/data/hapmap3.chr.7.fam\" \n \"/Users/huazhou/.julia/dev/OrdinalGWAS.jl/data/hapmap3.chr.8.bed\" \n \"/Users/huazhou/.julia/dev/OrdinalGWAS.jl/data/hapmap3.chr.8.bim\" \n \"/Users/huazhou/.julia/dev/OrdinalGWAS.jl/data/hapmap3.chr.8.fam\" \n \"/Users/huazhou/.julia/dev/OrdinalGWAS.jl/data/hapmap3.chr.9.bed\" \n \"/Users/huazhou/.julia/dev/OrdinalGWAS.jl/data/hapmap3.chr.9.bim\" \n \"/Users/huazhou/.julia/dev/OrdinalGWAS.jl/data/hapmap3.chr.9.fam\"Step 1: Fit the null model. Setting third argument plinkfile to nothing instructs ordinalgwas function to fit the null model only.nm = ordinalgwas(@formula(trait ~ sex), datadir * \"covariate.txt\", nothing)StatsModels.DataFrameRegressionModel{OrdinalMultinomialModel{Int64,Float64,LogitLink},Array{Float64,2}}\n\nFormula: trait ~ +sex\n\nCoefficients:\n      Estimate Std.Error  t value Pr(>|t|)\nθ1    -1.48564  0.358891 -4.13952    <1e-4\nθ2   -0.569479  0.341044 -1.66981   0.0959\nθ3    0.429815  0.339642  1.26549   0.2066\nβ1    0.424656  0.213914  1.98517   0.0480Step 2: GWAS for each chromosome.# this part can be submitted as separate jobs\nfor chr in 1:23\n    plinkfile = datadir * \"hapmap3.chr.\" * string(chr)\n    pvalfile = plinkfile * \".pval.txt\" \n    ordinalgwas(nm, plinkfile, pvalfile=pvalfile)\nend# show the result files\nreaddir(glob\"*.pval.txt\", datadir)23-element Array{String,1}:\n \"/Users/huazhou/.julia/dev/OrdinalGWAS.jl/data/hapmap3.chr.1.pval.txt\" \n \"/Users/huazhou/.julia/dev/OrdinalGWAS.jl/data/hapmap3.chr.10.pval.txt\"\n \"/Users/huazhou/.julia/dev/OrdinalGWAS.jl/data/hapmap3.chr.11.pval.txt\"\n \"/Users/huazhou/.julia/dev/OrdinalGWAS.jl/data/hapmap3.chr.12.pval.txt\"\n \"/Users/huazhou/.julia/dev/OrdinalGWAS.jl/data/hapmap3.chr.13.pval.txt\"\n \"/Users/huazhou/.julia/dev/OrdinalGWAS.jl/data/hapmap3.chr.14.pval.txt\"\n \"/Users/huazhou/.julia/dev/OrdinalGWAS.jl/data/hapmap3.chr.15.pval.txt\"\n \"/Users/huazhou/.julia/dev/OrdinalGWAS.jl/data/hapmap3.chr.16.pval.txt\"\n \"/Users/huazhou/.julia/dev/OrdinalGWAS.jl/data/hapmap3.chr.17.pval.txt\"\n \"/Users/huazhou/.julia/dev/OrdinalGWAS.jl/data/hapmap3.chr.18.pval.txt\"\n \"/Users/huazhou/.julia/dev/OrdinalGWAS.jl/data/hapmap3.chr.19.pval.txt\"\n \"/Users/huazhou/.julia/dev/OrdinalGWAS.jl/data/hapmap3.chr.2.pval.txt\" \n \"/Users/huazhou/.julia/dev/OrdinalGWAS.jl/data/hapmap3.chr.20.pval.txt\"\n \"/Users/huazhou/.julia/dev/OrdinalGWAS.jl/data/hapmap3.chr.21.pval.txt\"\n \"/Users/huazhou/.julia/dev/OrdinalGWAS.jl/data/hapmap3.chr.22.pval.txt\"\n \"/Users/huazhou/.julia/dev/OrdinalGWAS.jl/data/hapmap3.chr.23.pval.txt\"\n \"/Users/huazhou/.julia/dev/OrdinalGWAS.jl/data/hapmap3.chr.3.pval.txt\" \n \"/Users/huazhou/.julia/dev/OrdinalGWAS.jl/data/hapmap3.chr.4.pval.txt\" \n \"/Users/huazhou/.julia/dev/OrdinalGWAS.jl/data/hapmap3.chr.5.pval.txt\" \n \"/Users/huazhou/.julia/dev/OrdinalGWAS.jl/data/hapmap3.chr.6.pval.txt\" \n \"/Users/huazhou/.julia/dev/OrdinalGWAS.jl/data/hapmap3.chr.7.pval.txt\" \n \"/Users/huazhou/.julia/dev/OrdinalGWAS.jl/data/hapmap3.chr.8.pval.txt\" \n \"/Users/huazhou/.julia/dev/OrdinalGWAS.jl/data/hapmap3.chr.9.pval.txt\"In the rare situations where the multiple sets of Plink files lack the fam file or the corresponding bed and bim files have different filenames, users can explicitly supply bed filename, bim file name, and number of individuals. Replace Step 2 by Step 2\': GWAS for each chromosome.# this part can be submitted as separate jobs\nfor chr in 1:23\n    bedfile = datadir * \"hapmap3.chr.\" * string(chr) * \".bed\"\n    bimfile = datadir * \"hapmap3.chr.\" * string(chr) * \".bim\"\n    pvalfile = datadir * \"hapmap3.chr.\" * string(chr) * \".pval.txt\"\n    ordinalgwas(nm, bedfile, bimfile, 324; pvalfile=pvalfile)\nend# clean up\nisfile(\"ordinalgwas.null.txt\") && rm(\"ordinalgwas.null.txt\")\nisfile(datadir * \"fittednullmodel.jld2\") && rm(datadir * \"fittednullmodel.jld2\")\nfor chr in 1:23\n    pvalfile = datadir * \"hapmap3.chr.\" * string(chr) * \".pval.txt\"\n    isfile(pvalfile) && rm(pvalfile)\nend\nfor chr in 1:26\n    plinkfile = datadir * \"hapmap3.chr.\" * string(chr)\n    isfile(plinkfile * \".bed\") && rm(plinkfile * \".bed\")\n    isfile(plinkfile * \".fam\") && rm(plinkfile * \".fam\")\n    isfile(plinkfile * \".bim\") && rm(plinkfile * \".bim\")\nend"
},

{
    "location": "#Multiple-Plink-file-sets-on-cluster-1",
    "page": "OrdinalGWAS.jl",
    "title": "Multiple Plink file sets on cluster",
    "category": "section",
    "text": "We provide two scripts that successfully run on UCLA\'s Hoffman2 cluster using Julia v1.0.1 and PBS job schedulaer (qsub).The first script cluster_preparedata.jl creates a demo data set in current folder. Runjulia cluster_preparedata.jlon head node.run(`cat cluster_preparedata.jl`);#!/usr/local/bin/julia\n#\n# This script prepares a data set in current folder. \n# For each of chromosome 1-23, there is a set gzipped Plink files:\n# hapmap3.chr.1.bed.gz, hapmap3.chr.1.bim.gz, hapmap3.chr.1.fam.gz\n# hapmap3.chr.2.bed.gz, hapmap3.chr.2.bim.gz, hapmap3.chr.2.fam.gz\n# ...\n# hapmap3.chr.23.bed.gz, hapmap3.chr.23.bim.gz, hapmap3.chr.23.fam.gz\n# There is also a csv file \"covariate.txt\" that contains trait and covariates.\n#\n\n# install and load Julia packages\nusing Pkg\nif haskey(Pkg.installed(), \"SnpArrays\")\n    Pkg.update(\"SnpArrays\")\nelse\n    Pkg.add(PackageSpec(url=\"https://github.com/OpenMendel/SnpArrays.jl.git\"))\nend\nif haskey(Pkg.installed(), \"OrdinalMultinomialModels\")\n    Pkg.update(\"OrdinalMultinomialModels\")\nelse\n    Pkg.add(PackageSpec(url=\"https://github.com/OpenMendel/OrdinalMultinomialModels.jl.git\"))\nend\nif haskey(Pkg.installed(), \"OrdinalGWAS\")\n    Pkg.update(\"OrdinalGWAS\")\nelse\n    Pkg.add(PackageSpec(url=\"https://github.com/OpenMendel/OrdinalGWAS.jl.git\"))\nend\nusing OrdinalMultinomialModels, OrdinalGWAS, SnpArrays\n\n# split hapmap3 data according to chromosome\ndatadir = normpath(joinpath(dirname(pathof(OrdinalGWAS)), \"../data/\"))\nSnpArrays.split_plink(datadir * \"hapmap3\", :chromosome; prefix = \"hapmap3.chr.\")\n# compresse Plink files for chromosome 1-23\nfor chr in 1:23\n    plinkfile = \"hapmap3.chr.\" * string(chr)\n    SnpArrays.compress_plink(plinkfile)\nend\n# delete uncompressed chromosome Plink files\nfor chr in 1:26\n    plinkfile = \"hapmap3.chr.\" * string(chr)\n    isfile(plinkfile * \".bed\") && rm(plinkfile * \".bed\")\n    isfile(plinkfile * \".bim\") && rm(plinkfile * \".bim\")\n    isfile(plinkfile * \".fam\") && rm(plinkfile * \".fam\")\nend\n# copy covariate.txt file\ncp(datadir * \"covariate.txt\", joinpath(pwd(), \"covariate.txt\"))The second script cluster_run.jl first fits the null model then submits a separate job for each chromosome. Runjulia cluster_run.jlon head node.run(`cat cluster_run.jl`);#!/usr/local/bin/julia\n#\n# This script demonstrates how to submit multiple OrdinalGWAS runs from multiple sets of\n# Plink files on UCLA Hoffman2 cluster. It assumes that a demo data is available by\n# running `julia cluster_preparedata.jl` at current folder.\n#\n\nusing OrdinalGWAS, Serialization\n\n# Step 1: fit null model and save result to file `fittednullmodel.jls`\nnm = ordinalgwas(@formula(trait ~ sex), \"covariate.txt\", nothing)\nopen(\"fittednullmodel.jls\", \"w\") do io\n    Serialization.serialize(io, nm)\nend\n\n# Step 2: GWAS for each chromosome\nfor chr in 1:23\n    println(\"submit job for chromosome=$chr\")\n    jcode = \"using OrdinalGWAS, Serialization;\n    nm = open(deserialize, \\\"fittednullmodel.jls\\\");\n    bedfile = \\\"hapmap3.chr.\\\" * string($chr) * \\\".bed.gz\\\";\n    bimfile = \\\"hapmap3.chr.\\\" * string($chr) * \\\".bim.gz\\\";\n    pvalfile = \\\"hapmap3.chr.\\\" * string($chr) * \\\".pval.txt\\\";\n    ordinalgwas(nm, bedfile, bimfile, 324; pvalfile=pvalfile);\"\n    # prepare sh file for qsub\n    open(\"tmp.sh\", \"w\") do io\n        println(io, \"#!/bin/bash\")\n        println(io, \"#\\$ -cwd\")\n        println(io, \"# error = Merged with joblog\")\n        println(io, \"#\\$ -o joblog.\\$JOB_ID\")\n        println(io, \"#\\$ -j y\")\n        println(io, \"#\\$ -l h_rt=0:30:00,h_data=2G\") # request runtime and memory\n        println(io, \"#\\$ -pe shared 2\") # request # shared-memory nodes\n        println(io, \"# Email address to notify\")\n        println(io, \"#\\$ -M \\$USER@mail\")\n        println(io, \"# Notify when\")\n        println(io, \"#\\$ -m a\")\n        println(io)\n        println(io, \"# load the job environment:\")\n        println(io, \". /u/local/Modules/default/init/modules.sh\")\n        println(io, \"module load julia/1.0.1\") # available Julia version\n        println(io)\n        println(io, \"# run julia code\")\n        println(io, \"julia -e \'$jcode\' > output.\\$JOB_ID 2>&1\")\n    end\n    # submit job\n    run(`qsub tmp.sh`)\nend"
},

]}