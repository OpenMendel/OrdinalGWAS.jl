var documenterSearchIndex = {"docs": [

{
    "location": "#",
    "page": "PolrGWAS.jl",
    "title": "PolrGWAS.jl",
    "category": "page",
    "text": ""
},

{
    "location": "#PolrGWAS.jl-1",
    "page": "PolrGWAS.jl",
    "title": "PolrGWAS.jl",
    "category": "section",
    "text": "PolrGWAS.jl is a Julia package for performing genome-wide association studies (GWAS) for ordered categorical phenotypes using proportional odds model or ordred Probit model. It is useful when the phenotype takes ordered discrete values, e.g., disease status (undiagnosed, pre-disease, mild, moderate, severe). The package name follows the function polr in R package MASS."
},

{
    "location": "#Installation-1",
    "page": "PolrGWAS.jl",
    "title": "Installation",
    "category": "section",
    "text": "This package requires Julia v0.7.0 or later and two other unregistered packages SnpArrays and PolrModels. The package has not yet been registered and must be installed using the repository location. Start julia and use the ] key to switch to the package manager REPL(v1.0) pkg> add https://github.com/OpenMendel/SnpArrays.jl.git#juliav0.7\n(v1.0) pkg> add https://github.com/OpenMendel/PolrModels.jl.git\n(v1.0) pkg> add https://github.com/OpenMendel/PolrGWAS.jl.git# machine information for this tutorial\nversioninfo()Julia Version 1.0.2\nCommit d789231e99 (2018-11-08 20:11 UTC)\nPlatform Info:\n  OS: macOS (x86_64-apple-darwin14.5.0)\n  CPU: Intel(R) Core(TM) i7-6920HQ CPU @ 2.90GHz\n  WORD_SIZE: 64\n  LIBM: libopenlibm\n  LLVM: libLLVM-6.0.0 (ORCJIT, skylake)\nEnvironment:\n  JULIA_EDITOR = code# for use in this tutorial\nusing BenchmarkTools, CSV, PolrGWAS, SnpArrays"
},

{
    "location": "#Example-data-set-1",
    "page": "PolrGWAS.jl",
    "title": "Example data set",
    "category": "section",
    "text": "data folder of the package contains an example data set.;ls -l ../datatotal 3664\n-rw-r--r--  1 huazhou  staff     6844 Nov 23 17:58 covariate.txt\n-rw-r--r--  1 huazhou  staff  1128171 Nov 23 17:58 hapmap3.bed\n-rw-r--r--  1 huazhou  staff   388672 Nov 23 17:58 hapmap3.bim\n-rw-r--r--  1 huazhou  staff     7136 Nov 23 17:58 hapmap3.fam\n-rw-r--r--  1 huazhou  staff   332960 Nov 23 17:58 hapmap3.map\n-rw-r--r--  1 huazhou  staff      773 Nov 23 17:58 simtrait.jlcovariate.txt is a comma separated value (CSV) file containing the sample information, covariates sex, and phenotype trait. trait is coded as integer values 1, 2, 3 or 4. It was simulated from the script simtrait.jl. ;head ../data/covariate.txtfamid,perid,faid,moid,sex,trait\n2431,NA19916,0,0,1,4\n2424,NA19835,0,0,2,4\n2469,NA20282,0,0,2,4\n2368,NA19703,0,0,1,3\n2425,NA19901,0,0,2,3\n2427,NA19908,0,0,1,4\n2430,NA19914,0,0,2,4\n2470,NA20287,0,0,2,1\n2436,NA19713,0,0,2,3hapmap3 is a set of Plink files that contain the genotype information of samples. "
},

{
    "location": "#PolrGWAS.polrgwas",
    "page": "PolrGWAS.jl",
    "title": "PolrGWAS.polrgwas",
    "category": "function",
    "text": "polrgwas(nullformula, covfile, plinkfile)\npolrgwas(nullformula, df, plinkfile)\n\nPositional arguments\n\nnullformula::Formula: formula for the null model.\ncovfile::AbstractString: covariate file with one header line. One column    should be the ordered categorical phenotype coded as integers starting from 1.\ndf::DataFrame: DataFrame containing response and regressors.\nplinkfile::AbstractString: Plink file name without the bed, fam, or bim    extensions. If plinkfile==nothing, only null model is fitted.\n\nKeyword arguments\n\noutfile::AbstractString: output file prefix; default is polrgwas. Two CSV output files   prefix.nullmodel.txt and prefix.scoretest.txt (or prefix.lrttest.txt) will be written.\ncovtype::Vector{DataType}: type information for covfile. This is useful   when CSV.read(covarfile) has parsing errors.  \ntestformula::Formula: formula for test unit. Default is @formula(trait ~ 0 + snp).\ntest::Symbol: :score (default) or :LRT.  \nlink::GLM.Link: LogitLink() (default), ProbitLink(), CauchitLink(),   or CloglogLink().\nsnpmodel: ADDITIVE_MODEL (default), DOMINANT_MODEL, or RECESSIVE_MODEL.\ncolinds::Union{Nothing,AbstractVector{<:Integer}}: SNP indices.\nrowinds::Union{Nothing,AbstractVector{<:Integer}}: sample indices for bed file.\nsolver: an optimization solver supported by MathProgBase. Default is    NLoptSolver(algorithm=:LD_SLSQP, maxeval=4000). Another common choice is    IpoptSolver(print_level=0).\nverbose::Bool: default is false.\n\n\n\n\n\n"
},

{
    "location": "#Basic-usage-1",
    "page": "PolrGWAS.jl",
    "title": "Basic usage",
    "category": "section",
    "text": "The following command performs GWAS using the proportional odds model.polrgwas(@formula(trait ~ sex), \"../data/covariate.txt\", \"../data/hapmap3\")For documentation of the polrgwas function, type ?polrgwas in Julia REPL.polrgwas"
},

{
    "location": "#Formula-for-null-model-1",
    "page": "PolrGWAS.jl",
    "title": "Formula for null model",
    "category": "section",
    "text": "The first argument specifies the null model without SNP effects, e.g., @formula(trait ~ sex)."
},

{
    "location": "#Input-files-1",
    "page": "PolrGWAS.jl",
    "title": "Input files",
    "category": "section",
    "text": "polrgwas expects two input files: one for responses plus covariates (second argument), the other the Plink files for genotypes (third argument).Covariates and phenotype are available in a csv file covariate.txt, which has one header line for variable names. Variable trait is the ordered categorical phenotypes coded as integers 1 to 4. We want to include variable sex as the covariate in GWAS.;head ../data/covariate.txtfamid,perid,faid,moid,sex,trait\n2431,NA19916,0,0,1,4\n2424,NA19835,0,0,2,4\n2469,NA20282,0,0,2,4\n2368,NA19703,0,0,1,3\n2425,NA19901,0,0,2,3\n2427,NA19908,0,0,1,4\n2430,NA19914,0,0,2,4\n2470,NA20287,0,0,2,1\n2436,NA19713,0,0,2,3Genotype data is available as binary Plink files.;ls -l ../data/hapmap3.bed ../data/hapmap3.bim ../data/hapmap3.fam-rw-r--r--  1 huazhou  staff  1128171 Nov 23 17:58 ../data/hapmap3.bed\n-rw-r--r--  1 huazhou  staff   388672 Nov 23 17:58 ../data/hapmap3.bim\n-rw-r--r--  1 huazhou  staff     7136 Nov 23 17:58 ../data/hapmap3.famThere are 324 samples at 13,928 SNPs.size(SnpArray(\"../data/hapmap3.bed\"))(324, 13928)"
},

{
    "location": "#Output-files-1",
    "page": "PolrGWAS.jl",
    "title": "Output files",
    "category": "section",
    "text": "polrgwas outputs two files: polrgwas.nullmodel.txt and polrgwas.scoretest.txt. The prefix polrgwas of output files can be changed by the outfile keyword, e.g.,polrgwas(@formula(trait ~ sex), \"../data/covariate.txt\", \"../data/hapmap3\", outfile=\"hapmap3\")polrgwas.nullmodel.txt lists the estimated null model (without SNPs).  ;cat polrgwas.nullmodel.txtStatsModels.DataFrameRegressionModel{PolrModel{Int64,Float64,LogitLink},Array{Float64,2}}\n\nFormula: trait ~ +sex\n\nCoefficients:\n      Estimate Std.Error  t value Pr(>|t|)\nθ1    -1.48564  0.358891 -4.13952    <1e-4\nθ2   -0.569479  0.341044 -1.66981   0.0959\nθ3    0.429815  0.339642  1.26549   0.2066\nβ1    0.424656  0.213914  1.98517   0.0480polrgwas.scoretest.txt tallies the SNPs and their pvalues. ;head polrgwas.scoretest.txtchr,pos,snpid,maf,pval\n1,554484,rs10458597,0.0,1.0\n1,758311,rs12562034,0.07763975155279501,0.004565312839540994\n1,967643,rs2710875,0.32407407407407407,3.1082838285548695e-5\n1,1168108,rs11260566,0.19158878504672894,1.2168672367668912e-5\n1,1375074,rs1312568,0.441358024691358,0.008206860046175221\n1,1588771,rs35154105,0.0,1.0\n1,1789051,rs16824508,0.00462962962962965,0.5111981332544\n1,1990452,rs2678939,0.4537037037037037,0.29972829571847825\n1,2194615,rs7553178,0.22685185185185186,0.1713331245805063"
},

{
    "location": "#Input-non-genetic-data-as-DataFrame-1",
    "page": "PolrGWAS.jl",
    "title": "Input non-genetic data as DataFrame",
    "category": "section",
    "text": "Internally polrgwas parses the covariate file as a DataFrame by CSV.read(covfile). For covariate file of other format, users can parse first and then input a DataFrame to polrgwas directly.polrgwas(@formula(trait ~ sex), df, \"../data/hapmap3\")note: Note\nUsers should always make sure that individuals in covariate file or DataFrame match those in Plink fam file. For example, following code checks that the first 2 columns of the covariate.txt file match the first 2 columns of the hapmap3.fam file exactly.covdf = CSV.read(\"../data/covariate.txt\")\nplkfam = CSV.read(\"../data/hapmap3.fam\", header=0, delim=\' \')\nall(covdf[1] .== plkfam[1]) && all(covdf[2] .== plkfam[2])true"
},

{
    "location": "#Timing-1",
    "page": "PolrGWAS.jl",
    "title": "Timing",
    "category": "section",
    "text": "For this moderate-sized data set, polrgwas takes less than 0.2 second.@btime(polrgwas(@formula(trait ~ sex), \"../data/covariate.txt\", \"../data/hapmap3\"))  129.735 ms (639113 allocations: 29.12 MiB)# clean up\nrm(\"polrgwas.scoretest.txt\")\nrm(\"polrgwas.nullmodel.txt\")"
},

{
    "location": "#Link-functions-1",
    "page": "PolrGWAS.jl",
    "title": "Link functions",
    "category": "section",
    "text": "The link keyword argument of polrgwas can take value:  LogitLink(), proportional odds model (default),  \nProbitLink(), ordred Probit model,  \nCloglogLink(), proportional hazards model), or \nCauchyLink().For example, to perform GWAS using the ordred Probit modelpolrgwas(@formula(trait ~ sex), \"../data/covariate.txt\", \"../data/hapmap3\", \n    link=ProbitLink(), outfile=\"opm\")The estimates in null model and p-values are slightly different from proportional odds moodel.;cat opm.nullmodel.txtStatsModels.DataFrameRegressionModel{PolrModel{Int64,Float64,ProbitLink},Array{Float64,2}}\n\nFormula: trait ~ +sex\n\nCoefficients:\n      Estimate Std.Error  t value Pr(>|t|)\nθ1   -0.866156  0.210677 -4.11129    <1e-4\nθ2   -0.359878  0.205817 -1.74854   0.0813\nθ3    0.247054  0.205382   1.2029   0.2299\nβ1    0.251058  0.128225  1.95795   0.0511;head opm.scoretest.txtchr,pos,snpid,maf,pval\n1,554484,rs10458597,0.0,1.0\n1,758311,rs12562034,0.07763975155279501,0.010076916742300138\n1,967643,rs2710875,0.32407407407407407,2.6272564941853933e-5\n1,1168108,rs11260566,0.19158878504672894,1.0897484851078458e-5\n1,1375074,rs1312568,0.441358024691358,0.005102883990438149\n1,1588771,rs35154105,0.0,1.0\n1,1789051,rs16824508,0.00462962962962965,0.48653776297859236\n1,1990452,rs2678939,0.4537037037037037,0.33231290090455434\n1,2194615,rs7553178,0.22685185185185186,0.25915513977197435rm(\"opm.nullmodel.txt\")\nrm(\"opm.scoretest.txt\")"
},

{
    "location": "#SNP-models-1",
    "page": "PolrGWAS.jl",
    "title": "SNP models",
    "category": "section",
    "text": "Genotypes are translated into numeric values according to different genetic model, which is specified by the snpmodel keyword. Default is ADDITIVE_MODEL.Genotype SnpArray ADDITIVE_MODEL DOMINANT_MODEL RECESSIVE_MODEL\nA1,A1 0x00 0 0 0\nmissing 0x01 NaN NaN NaN\nA1,A2 0x02 1 1 0\nA2,A2 0x03 2 1 1"
},

{
    "location": "#SNP-and/or-sample-masks-1",
    "page": "PolrGWAS.jl",
    "title": "SNP and/or sample masks",
    "category": "section",
    "text": "In practice, we often perform GWAS on selected SNPs and/or selected samples. They can be specified by the colinds and rowinds keywords of polrgwas function.For example, to perform GWAS on SNPs with minor allele frequency (MAF) above 0.05# create SNP mask\nsnpinds = maf(SnpArray(\"../data/hapmap3.bed\")) .≥ 0.0513928-element BitArray{1}:\n false\n  true\n  true\n  true\n  true\n false\n false\n  true\n  true\n  true\n false\n  true\n  true\n     ⋮\n  true\n  true\n  true\n  true\n  true\n  true\n  true\n  true\n  true\n false\n false\n false# GWAS on selected SNPs\n@time polrgwas(@formula(trait ~ sex), \"../data/covariate.txt\", \"../data/hapmap3\", \n    colinds = snpinds, outfile=\"commonvariant\")  0.222907 seconds (767.25 k allocations: 36.103 MiB, 6.67% gc time);cat commonvariant.nullmodel.txtStatsModels.DataFrameRegressionModel{PolrModel{Int64,Float64,LogitLink},Array{Float64,2}}\n\nFormula: trait ~ +sex\n\nCoefficients:\n      Estimate Std.Error  t value Pr(>|t|)\nθ1    -1.48564  0.358891 -4.13952    <1e-4\nθ2   -0.569479  0.341044 -1.66981   0.0959\nθ3    0.429815  0.339642  1.26549   0.2066\nβ1    0.424656  0.213914  1.98517   0.0480;head commonvariant.scoretest.txtchr,pos,snpid,maf,pval\n1,758311,rs12562034,0.07763975155279501,0.004565312839540994\n1,967643,rs2710875,0.32407407407407407,3.1082838285548695e-5\n1,1168108,rs11260566,0.19158878504672894,1.2168672367668912e-5\n1,1375074,rs1312568,0.441358024691358,0.008206860046175221\n1,1990452,rs2678939,0.4537037037037037,0.29972829571847825\n1,2194615,rs7553178,0.22685185185185186,0.1713331245805063\n1,2396747,rs13376356,0.1448598130841121,0.5320416198875456\n1,2823603,rs1563468,0.4830246913580247,0.225191391783573\n1,3025087,rs6690373,0.2538699690402477,0.7018469417717486# extra headline in commonvariant.scoretest.txt\ncountlines(\"commonvariant.scoretest.txt\"), count(snpinds)(12086, 12085)# clean up\nrm(\"commonvariant.scoretest.txt\")\nrm(\"commonvariant.nullmodel.txt\")User should be particularly careful when using the rowinds keyword. Selected rows in SnpArray should exactly match the samples in the null model. Otherwise the results are meaningless."
},

{
    "location": "#Likelihood-ratio-test-(LRT)-1",
    "page": "PolrGWAS.jl",
    "title": "Likelihood ratio test (LRT)",
    "category": "section",
    "text": "By default, polrgwas calculates p-value for each SNP using score test. Score test is fast because it doesn\'t require fitting alternative model for each SNP. User can request likelihood ratio test (LRT) using keyword test=:LRT. LRT is much slower but may be more powerful than score test.@time polrgwas(@formula(trait ~ sex), \"../data/covariate.txt\", \"../data/hapmap3\", \n    test=:LRT, outfile=\"lrt\") 19.855754 seconds (8.19 M allocations: 2.045 GiB, 2.01% gc time)Test result is output to outfile.lrttest.txt file;head lrt.lrttest.txtchr,pos,snpid,maf,effect,pval\n1,554484,rs10458597,0.0,0.0,1.0\n1,758311,rs12562034,0.07763975155279501,-1.0057833719544331,0.0019185836579804134\n1,967643,rs2710875,0.32407407407407407,-0.6488560566295055,1.805050556976241e-5\n1,1168108,rs11260566,0.19158878504672894,-0.9157225669357879,5.873384712685568e-6\n1,1375074,rs1312568,0.441358024691358,-0.33181366525772593,0.008081022577832324\n1,1588771,rs35154105,0.0,0.0,1.0\n1,1789051,rs16824508,0.00462962962962965,-0.7338026388701573,0.5169027130129711\n1,1990452,rs2678939,0.4537037037037037,-0.13586499231819726,0.29946402200912603\n1,2194615,rs7553178,0.22685185185185186,-0.2512075640440123,0.16151069094439868Note the extra effect column, which is the effect size (regression coefficient) for each SNP. # clean up\nrm(\"lrt.lrttest.txt\")\nrm(\"lrt.nullmodel.txt\")In this example, GWAS by score test takes less than 0.2 second, while GWAS by LRT takes about 20 seconds. About 100 fold difference in run time. "
},

{
    "location": "#Score-test-for-screening,-LRT-for-power-1",
    "page": "PolrGWAS.jl",
    "title": "Score test for screening, LRT for power",
    "category": "section",
    "text": "For large data sets, a practical solution is to perform score test first, then re-do LRT for the most promising SNPs according to score test p-values.Step 1: Perform score test GWAS, results in hapmap.scoretest.txt.@time polrgwas(@formula(trait ~ sex), \"../data/covariate.txt\", \"../data/hapmap3\", \n    test=:score, outfile=\"hapmap\", verbose=false)  0.232228 seconds (717.17 k allocations: 33.237 MiB, 6.99% gc time);head hapmap.scoretest.txtchr,pos,snpid,maf,pval\n1,554484,rs10458597,0.0,1.0\n1,758311,rs12562034,0.07763975155279501,0.004565312839540994\n1,967643,rs2710875,0.32407407407407407,3.1082838285548695e-5\n1,1168108,rs11260566,0.19158878504672894,1.2168672367668912e-5\n1,1375074,rs1312568,0.441358024691358,0.008206860046175221\n1,1588771,rs35154105,0.0,1.0\n1,1789051,rs16824508,0.00462962962962965,0.5111981332544\n1,1990452,rs2678939,0.4537037037037037,0.29972829571847825\n1,2194615,rs7553178,0.22685185185185186,0.1713331245805063Step 2: Sort score test p-values and find top 10 SNPs.scorepvals = CSV.read(\"hapmap.scoretest.txt\")[5]13928-element Array{Union{Missing, Float64},1}:\n 1.0                  \n 0.004565312839540994 \n 3.1082838285548695e-5\n 1.2168672367668912e-5\n 0.008206860046175221 \n 1.0                  \n 0.5111981332544      \n 0.29972829571847825  \n 0.1713331245805063   \n 0.5320416198875456   \n 1.0                  \n 0.225191391783573    \n 0.7018469417717486   \n ⋮                    \n 0.3659920717810955   \n 0.5612543378640146   \n 0.11646522388490779  \n 0.5216852471858483   \n 0.44711860378438606  \n 0.16565527326637602  \n 0.65590930553447     \n 0.89481895007441     \n 0.6559093055344556   \n 0.3462815220311233   \n 0.2316733923395225   \n 0.2651082379376032tophits = sortperm(scorepvals)[1:10]10-element Array{Int64,1}:\n  6063\n  5071\n  2458\n     4\n  3291\n 13737\n  5293\n     3\n  6256\n  4183scorepvals[tophits]10-element Array{Union{Missing, Float64},1}:\n 1.3080149099181335e-6\n 6.536722765052079e-6 \n 9.664742185669054e-6 \n 1.2168672367668912e-5\n 1.802746001833127e-5 \n 2.0989542284213636e-5\n 2.6844521269963608e-5\n 3.1082838285548695e-5\n 4.1010912875160476e-5\n 4.2966265138454806e-5Step 3: Re-do LRT on top hits.@time polrgwas(@formula(trait ~ sex), \"../data/covariate.txt\", \"../data/hapmap3\", \n    colinds=tophits, test=:LRT, outfile=\"hapmap\")  0.159806 seconds (351.95 k allocations: 18.710 MiB, 4.96% gc time);head hapmap.lrttest.txtchr,pos,snpid,maf,effect,pval\n1,967643,rs2710875,0.32407407407407407,-0.6488560566295055,1.805050556976241e-5\n1,1168108,rs11260566,0.19158878504672894,-0.9157225669357879,5.873384712685568e-6\n3,36821790,rs4678553,0.23456790123456794,0.7424952268973518,1.1303825016262592e-5\n4,11017683,rs16881446,0.27554179566563464,-0.7870581482955515,1.1105427468799613e-5\n5,3739190,rs12521166,0.0679012345679012,1.1468852997925316,4.781288229657399e-5\n6,7574576,rs1885466,0.17746913580246915,0.8750621092263019,7.272346896740631e-6\n6,52474721,rs2073183,0.1826625386996904,0.7790794914858663,5.069394513906121e-5\n7,41152376,rs28880,0.3379629629629629,-0.814633902445351,9.180126530294943e-7\n7,84223996,rs4128623,0.07870370370370372,1.0022229316338573,6.587895464657512e-5# clean up\nrm(\"hapmap.nullmodel.txt\")\nrm(\"hapmap.lrttest.txt\")\nrm(\"hapmap.scoretest.txt\")"
},

{
    "location": "#GxE-or-other-interactions-1",
    "page": "PolrGWAS.jl",
    "title": "GxE or other interactions",
    "category": "section",
    "text": "In many applications, we want to test SNP effect and/or its interaction with other terms. testformula keyword specifies the test unit besides the covariates in nullformula. In following example, keyword testformula=@formula(trait ~ snp + snp & sex) instructs polrgwas to test joint effect of snp and snp & sex interaction.polrgwas(@formula(trait ~ sex), \"../data/covariate.txt\", \"../data/hapmap3\", \n    outfile=\"GxE\", testformula=@formula(trait ~ snp + snp & sex));head GxE.scoretest.txtchr,pos,snpid,maf,pval\n1,554484,rs10458597,0.0,1.0\n1,758311,rs12562034,0.07763975155279501,0.017446010412254197\n1,967643,rs2710875,0.32407407407407407,0.0001667073239489097\n1,1168108,rs11260566,0.19158878504672894,4.763762457893366e-5\n1,1375074,rs1312568,0.441358024691358,0.029138471242993652\n1,1588771,rs35154105,0.0,1.0\n1,1789051,rs16824508,0.00462962962962965,0.2964363114944328\n1,1990452,rs2678939,0.4537037037037037,0.37924580479348785\n1,2194615,rs7553178,0.22685185185185186,0.325582269932396# clean up\nrm(\"GxE.nullmodel.txt\")\nrm(\"GxE.scoretest.txt\")To test the linear and quadratic SNP effects jointlypolrgwas(@formula(trait ~ 0 + sex), \"../data/covariate.txt\", \"../data/hapmap3\", \n    outfile=\"quadratic\", testformula=@formula(trait ~ snp + snp & snp));head quadratic.scoretest.txtchr,pos,snpid,maf,pval\n1,554484,rs10458597,0.0,1.0\n1,758311,rs12562034,0.07763975155279501,1.0\n1,967643,rs2710875,0.32407407407407407,1.0\n1,1168108,rs11260566,0.19158878504672894,1.0\n1,1375074,rs1312568,0.441358024691358,1.0\n1,1588771,rs35154105,0.0,1.0\n1,1789051,rs16824508,0.00462962962962965,1.0\n1,1990452,rs2678939,0.4537037037037037,1.0\n1,2194615,rs7553178,0.22685185185185186,1.0# clean up\nrm(\"quadratic.nullmodel.txt\")\nrm(\"quadratic.scoretest.txt\")"
},

{
    "location": "#Plots-1",
    "page": "PolrGWAS.jl",
    "title": "Plots",
    "category": "section",
    "text": "To plot the GWAS results, use the MendelPlots (TODO) package."
},

{
    "location": "#Docker-1",
    "page": "PolrGWAS.jl",
    "title": "Docker",
    "category": "section",
    "text": "For ease of using PolrGWAS, we provide a Dockerfile so users don\'t need to install Julia and required packages. Only Docker app needs to be installed in order to run analysis. Following is tested on Docker 2.0.0.0-mac78.Step 1: Create a Dockerfile with content shown below. You can copy from here.;cat ../docker/DockerfileFROM centos:7\n\nWORKDIR /root\n\nRUN yum update -y && yum install -y epel-release && yum clean all\n\nRUN yum update -y && yum install -y \\\n    cmake \\\n    curl-devel \\\n    expat-devel \\\n    gcc \\\n    gcc-c++ \\\n    gcc-gfortran \\\n    gettext-devel \\\n    make \\\n    openssl \\\n    openssl098e \\\n    openssl-devel \\\n    patch \\\n    svn \\\n    wget \\\n    which \\\n    && yum clean all    \n\nENV PATH /usr/local/sbin:/usr/local/bin:$PATH\n\nENV LD_LIBRARY_PATH /usr/local/lib:/usr/local/lib64\n\n# GIT - https://git-scm.com/\n# http://tecadmin.net/install-git-2-0-on-centos-rhel-fedora/#\nENV GIT_VER 2.19.1\n\nRUN wget https://www.kernel.org/pub/software/scm/git/git-$GIT_VER.tar.gz \\\n    && tar xf git-$GIT_VER.tar.gz && cd git-$GIT_VER \\\n    && make -j\"$(nproc --all)\" prefix=/usr/local all \\\n    && make prefix=/usr/local -j\"$(nproc --all)\" install \\\n    && cd .. && rm -f git-$GIT_VER.tar.gz && rm -rf git-$GIT_VER\n\n# Makes git use https by default\nRUN git config --global url.\"https://\".insteadOf git://\n\n# Julia\nENV JULIA_VER_MAJ 1.0\nENV JULIA_VER_MIN .1\nENV JULIA_VER $JULIA_VER_MAJ$JULIA_VER_MIN\n\nRUN wget https://julialang-s3.julialang.org/bin/linux/x64/$JULIA_VER_MAJ/julia-$JULIA_VER-linux-x86_64.tar.gz \\\n    && mkdir /usr/local/julia \\\n    && tar xf julia-$JULIA_VER-linux-x86_64.tar.gz --directory /usr/local/julia --strip-components=1 \\\n    && ln -s /usr/local/julia/bin/julia /usr/local/bin/julia \\\n    && rm -f julia-$JULIA_VER-linux-x86_64.tar.gz\n\nENV JULIA_PKGDIR /usr/local/julia/share/julia/site\n\nRUN julia -e \'using Pkg; \\\n    Pkg.add([ \\\n    PackageSpec(url=\"https://github.com/OpenMendel/SnpArrays.jl.git\", rev=\"juliav0.7\"), \\\n    PackageSpec(url=\"https://github.com/OpenMendel/PolrModels.jl.git\", rev=\"master\"), \\\n    PackageSpec(url=\"https://github.com/OpenMendel/PolrGWAS.jl.git\", rev=\"master\") \\\n    ]); \\\n    Pkg.test(\"PolrGWAS\");\'Step 2: Build a docker image polrgwas-app, assuming that the Dockerfile is located at ../docker folder. Building the image for the first time can take up to 10 minutes; but it only needs to be done once.;docker build -t polrgwas-app ../docker/Sending build context to Docker daemon  3.584kB\nStep 1/15 : FROM centos:7\n ---> 75835a67d134\nStep 2/15 : WORKDIR /root\n ---> Using cache\n ---> f052f835f42d\nStep 3/15 : RUN yum update -y && yum install -y epel-release && yum clean all\n ---> Using cache\n ---> a75b6baa5581\nStep 4/15 : RUN yum update -y && yum install -y     cmake     curl-devel     expat-devel     gcc     gcc-c++     gcc-gfortran     gettext-devel     make     openssl     openssl098e     openssl-devel     patch     svn     wget     which     && yum clean all\n ---> Using cache\n ---> dcf0f6ba6334\nStep 5/15 : ENV PATH /usr/local/sbin:/usr/local/bin:$PATH\n ---> Using cache\n ---> 7861f46c17eb\nStep 6/15 : ENV LD_LIBRARY_PATH /usr/local/lib:/usr/local/lib64\n ---> Using cache\n ---> 4d3222181efc\nStep 7/15 : ENV GIT_VER 2.19.1\n ---> Using cache\n ---> b23c5cb67786\nStep 8/15 : RUN wget https://www.kernel.org/pub/software/scm/git/git-$GIT_VER.tar.gz     && tar xf git-$GIT_VER.tar.gz && cd git-$GIT_VER     && make -j\"$(nproc --all)\" prefix=/usr/local all     && make prefix=/usr/local -j\"$(nproc --all)\" install     && cd .. && rm -f git-$GIT_VER.tar.gz && rm -rf git-$GIT_VER\n ---> Using cache\n ---> 6beaafd797ec\nStep 9/15 : RUN git config --global url.\"https://\".insteadOf git://\n ---> Using cache\n ---> f1aec5dd453f\nStep 10/15 : ENV JULIA_VER_MAJ 1.0\n ---> Using cache\n ---> 42eb17d93684\nStep 11/15 : ENV JULIA_VER_MIN .1\n ---> Using cache\n ---> 4289e3f44a22\nStep 12/15 : ENV JULIA_VER $JULIA_VER_MAJ$JULIA_VER_MIN\n ---> Using cache\n ---> 2e9b9213c88c\nStep 13/15 : RUN wget https://julialang-s3.julialang.org/bin/linux/x64/$JULIA_VER_MAJ/julia-$JULIA_VER-linux-x86_64.tar.gz     && mkdir /usr/local/julia     && tar xf julia-$JULIA_VER-linux-x86_64.tar.gz --directory /usr/local/julia --strip-components=1     && ln -s /usr/local/julia/bin/julia /usr/local/bin/julia     && rm -f julia-$JULIA_VER-linux-x86_64.tar.gz\n ---> Using cache\n ---> 2510e817738a\nStep 14/15 : ENV JULIA_PKGDIR /usr/local/julia/share/julia/site\n ---> Using cache\n ---> 3692e40f50fa\nStep 15/15 : RUN julia -e \'using Pkg;     Pkg.add([     PackageSpec(url=\"https://github.com/OpenMendel/SnpArrays.jl.git\", rev=\"juliav0.7\"),     PackageSpec(url=\"https://github.com/OpenMendel/PolrModels.jl.git\", rev=\"master\"),     PackageSpec(url=\"https://github.com/OpenMendel/PolrGWAS.jl.git\", rev=\"master\")     ]);     Pkg.test(\"PolrGWAS\");\'\n ---> Using cache\n ---> af7e8765610c\nSuccessfully built af7e8765610c\nSuccessfully tagged polrgwas-app:latestStep 3: Suppose data files are located at /Users/huazhou/.julia/dev/PolrGWAS/data folder, run analysis by;docker run -v /Users/huazhou/.julia/dev/PolrGWAS/data:/data -t polrgwas-app julia -e \'using PolrGWAS; polrgwas(@formula(trait ~ sex), \"/data/covariate.txt\", \"/data/hapmap3\", outfile=\"/data/polrgwas\");\'Here  -t polrgwas-app creates a container using the polrgwas-app image build in step 2.  \n-v /Users/huazhou/.julia/dev/PolrGWAS/data:/data maps the /Users/huazhou/.julia/dev/PolrGWAS/data folder on host machine to the /data folder within the container. \njulia -e \'using PolrGWAS; polrgwas(@formula(trait ~ 0 + sex), \"/data/covariate.txt\", \"/data/hapmap3\", outfile=\"/data/polrgwas\"); calls Julia and runs polrgwas function. The output files are at /Users/huazhou/.julia/dev/PolrGWAS/data. ;ls -l /Users/huazhou/.julia/dev/PolrGWAS/data/total 5320\n-rw-r--r--  1 huazhou  staff     6844 Nov 23 17:58 covariate.txt\n-rw-r--r--  1 huazhou  staff  1128171 Nov 23 17:58 hapmap3.bed\n-rw-r--r--  1 huazhou  staff   388672 Nov 23 17:58 hapmap3.bim\n-rw-r--r--  1 huazhou  staff     7136 Nov 23 17:58 hapmap3.fam\n-rw-r--r--  1 huazhou  staff   332960 Nov 23 17:58 hapmap3.map\n-rw-r--r--  1 huazhou  staff      347 Dec  1 23:06 polrgwas.nullmodel.txt\n-rw-r--r--  1 huazhou  staff   842401 Dec  1 23:06 polrgwas.scoretest.txt\n-rw-r--r--  1 huazhou  staff      773 Nov 23 17:58 simtrait.jl# clean up\nrm(\"/Users/huazhou/.julia/dev/PolrGWAS/data/polrgwas.nullmodel.txt\")\nrm(\"/Users/huazhou/.julia/dev/PolrGWAS/data/polrgwas.scoretest.txt\")"
},

]}
