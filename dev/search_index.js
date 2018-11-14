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
    "text": "PolrGWAS.jl is a Julia package for performing genome-wide association studies (GWAS) for ordered categorical phenotypes. It is useful when the phenotype takes ordered discrete values, e.g., disease status (undiagnosed, pre-disease, mild, moderate, severe)."
},

{
    "location": "#Installation-1",
    "page": "PolrGWAS.jl",
    "title": "Installation",
    "category": "section",
    "text": "This package requires Julia v0.7.0 or later and two other unregistered packages SnpArrays and PolrModels. The package has not yet been registered and must be installed using the repository location. Start julia and use the ] key to switch to the package manager REPL(v0.7) pkg> add https://github.com/OpenMendel/SnpArrays.git#juliav0.7\n(v0.7) pkg> add https://github.com/OpenMendel/PolrModels.git\n(v0.7) pkg> add https://github.com/OpenMendel/PolrGWAS.gitversioninfo()Julia Version 0.7.0\nCommit a4cb80f3ed (2018-08-08 06:46 UTC)\nPlatform Info:\n  OS: macOS (x86_64-apple-darwin14.5.0)\n  CPU: Intel(R) Core(TM) i7-6920HQ CPU @ 2.90GHz\n  WORD_SIZE: 64\n  LIBM: libopenlibm\n  LLVM: libLLVM-6.0.0 (ORCJIT, skylake)\nEnvironment:\n  JULIA_EDITOR = code# for use in this tutorial\nusing BenchmarkTools, CSV, PolrGWAS, SnpArrays"
},

{
    "location": "#Example-data-set-1",
    "page": "PolrGWAS.jl",
    "title": "Example data set",
    "category": "section",
    "text": "data folder of the package contains an example data set.;ls -l ../datatotal 3664\n-rw-r--r--  1 huazhou  staff     6844 Nov 11 21:09 covariate.txt\n-rw-r--r--  1 huazhou  staff  1128171 Nov 11 21:09 hapmap3.bed\n-rw-r--r--  1 huazhou  staff   388672 Nov 11 21:09 hapmap3.bim\n-rw-r--r--  1 huazhou  staff     7136 Nov 11 21:09 hapmap3.fam\n-rw-r--r--  1 huazhou  staff   332960 Nov 11 21:09 hapmap3.map\n-rw-r--r--  1 huazhou  staff      773 Nov 11 21:09 simtrait.jlcovariate.txt is a comma separated value (CSV) file containing the sample information, covariates sex, and phenotype trait. trait is coded as integer values 1, 2, 3 or 4. It was simulated from the script simtrait.jl. ;head -20 ../data/covariate.txtfamid,perid,faid,moid,sex,trait\n2431,NA19916,0,0,1,4\n2424,NA19835,0,0,2,4\n2469,NA20282,0,0,2,4\n2368,NA19703,0,0,1,3\n2425,NA19901,0,0,2,3\n2427,NA19908,0,0,1,4\n2430,NA19914,0,0,2,4\n2470,NA20287,0,0,2,1\n2436,NA19713,0,0,2,3\n2426,NA19904,0,0,1,1\n2431,NA19917,0,0,2,1\n2436,NA19982,0,0,1,2\n2487,NA20340,0,0,1,4\n2427,NA19909,0,0,2,4\n2424,NA19834,0,0,1,4\n2480,NA20317,0,0,2,4\n2418,NA19818,0,0,1,1\n2490,NA20346,0,0,1,2\n2433,NA19921,0,0,2,4hapmap3 is a set of Plink files that contain the genotype information of samples. "
},

{
    "location": "#PolrGWAS.polrgwas",
    "page": "PolrGWAS.jl",
    "title": "PolrGWAS.polrgwas",
    "category": "function",
    "text": "polrgwas(nullformula, covfile, plkfile)\npolrgwas(nullformula, df, plkfile)\n\nPositional arguments\n\nnullformula::Formula: formula for the null model.\ncovfile::AbstractString: covariate file with one header line. One column    should be the ordered categorical phenotype coded as integers starting from 1.\ndf::DataFrame: DataFrame containing response and regressors.\nplkfile::AbstractString: Plink file name without the bed, fam, or bim    extensions. If plkfile==nothing, only null model is fitted.\n\nKeyword arguments\n\noutfile::AbstractString: output file prefix; default is polrgwas. Two CSV output files   prefix.nullmodel.txt and prefix.scoretest.txt (or prefix.lrttest.txt) will be written.\ncovtype::Vector{DataType}: type information for covfile. This is useful   when CSV.read(covarfile) has parsing errors.  \ntestformula::Formula: formula for test unit. Default is @formula(~ 0 + snp).\ntest::Symbol: :score (default) or :LRT.\nlink::GLM.Link: LogitLink() (default), ProbitLink(), CauchitLink(),   or CloglogLink().\ncolinds::Union{Nothing,AbstractVector{<:Integer}}: SNP indices.\nrowinds::Union{Nothing,AbstractVector{<:Integer}}: sample indices for bed file.\nsolver: an optimization solver supported by MathProgBase. Default is    NLoptSolver(algorithm=:LD_SLSQP, maxeval=4000). Another common choice is    IpoptSolver(print_level=0).\nverbose::Bool: default is false.\n\n\n\n\n\n"
},

{
    "location": "#Basic-usage-1",
    "page": "PolrGWAS.jl",
    "title": "Basic usage",
    "category": "section",
    "text": "The following command performs GWAS using the proportional odds model.polrgwas(@formula(trait ~ 0 + sex), \"../data/covariate.txt\", \"../data/hapmap3\")For documentation of the polrgwas function, type ?polrgwas in Julia REPL.polrgwas"
},

{
    "location": "#Formula-for-null-model-1",
    "page": "PolrGWAS.jl",
    "title": "Formula for null model",
    "category": "section",
    "text": "The first argument specifies the null model without SNP effects, e.g., @formula(trait ~ 0 + sex). It is important to exclude the intercept because proportional odds model automatically incorporates intercepts for modeling purpose."
},

{
    "location": "#Input-files-1",
    "page": "PolrGWAS.jl",
    "title": "Input files",
    "category": "section",
    "text": "polrgwas expects two input files: one for responses and covariates (second argument), the other the Plink files for genotypes (third argument).Covariates and phenotype are available in a csv file covariate.txt, which has one header line for variable names. Variable trait is the ordered categorical phenotypes coded as integers 1 to 4. We want to include variable sex as the covariate in GWAS.;head -20 ../data/covariate.txtfamid,perid,faid,moid,sex,trait\n2431,NA19916,0,0,1,4\n2424,NA19835,0,0,2,4\n2469,NA20282,0,0,2,4\n2368,NA19703,0,0,1,3\n2425,NA19901,0,0,2,3\n2427,NA19908,0,0,1,4\n2430,NA19914,0,0,2,4\n2470,NA20287,0,0,2,1\n2436,NA19713,0,0,2,3\n2426,NA19904,0,0,1,1\n2431,NA19917,0,0,2,1\n2436,NA19982,0,0,1,2\n2487,NA20340,0,0,1,4\n2427,NA19909,0,0,2,4\n2424,NA19834,0,0,1,4\n2480,NA20317,0,0,2,4\n2418,NA19818,0,0,1,1\n2490,NA20346,0,0,1,2\n2433,NA19921,0,0,2,4Genotype data is available as binary Plink files.;ls -l ../data/hapmap3.bed ../data/hapmap3.bim ../data/hapmap3.fam-rw-r--r--  1 huazhou  staff  1128171 Nov 11 21:09 ../data/hapmap3.bed\n-rw-r--r--  1 huazhou  staff   388672 Nov 11 21:09 ../data/hapmap3.bim\n-rw-r--r--  1 huazhou  staff     7136 Nov 11 21:09 ../data/hapmap3.famThere are 324 samples at 13,928 SNPs.size(SnpArray(\"../data/hapmap3.bed\"))(324, 13928)"
},

{
    "location": "#Output-files-1",
    "page": "PolrGWAS.jl",
    "title": "Output files",
    "category": "section",
    "text": "polrgwas outputs two files: polrgwas.nullmodel.txt and polrgwas.scoretest.txt. polrgwas.nullmodel.txt lists the estimated null model (without SNPs).  ;cat polrgwas.nullmodel.txtStatsModels.DataFrameRegressionModel{PolrModel{Int64,Float64,LogitLink},Array{Float64,2}}\n\nFormula: trait ~ +sex\n\nCoefficients:\n      Estimate Std.Error  t value Pr(>|t|)\nθ1    -1.48564  0.358713 -4.14157    <1e-4\nθ2   -0.569479  0.340649 -1.67175   0.0956\nθ3    0.429815  0.339266   1.2669   0.2061\nβ1    0.424656  0.213911   1.9852   0.0480polrgwas.scoretest.txt tallies the SNPs and their pvalues. ;head polrgwas.scoretest.txtchr,pos,snpid,maf,pval\n1,554484,rs10458597,0.0,1.0\n1,758311,rs12562034,0.07763975155279501,0.003001230754791795\n1,967643,rs2710875,0.32407407407407407,2.5117214960313727e-5\n1,1168108,rs11260566,0.19158878504672894,1.13731122530895e-5\n1,1375074,rs1312568,0.441358024691358,0.008317358366815325\n1,1588771,rs35154105,0.0,1.0\n1,1789051,rs16824508,0.00462962962962965,0.5274428530031774\n1,1990452,rs2678939,0.4537037037037037,0.29988429741739986\n1,2194615,rs7553178,0.22685185185185186,0.16436415589171793The prefix of output files can be changed by the outfile keyword, e.g.,polrgwas(@formula(trait ~ 0 + sex), \"../data/covariate.txt\", \"../data/hapmap3\", \n    outfile=\"hapmap3\")"
},

{
    "location": "#Input-non-genetic-data-as-DataFrame-1",
    "page": "PolrGWAS.jl",
    "title": "Input non-genetic data as DataFrame",
    "category": "section",
    "text": "Internally polrgwas parses the covariate file as a DataFrame by CSV.read(covfile). For covariate file of other format, users can parse first and then input a DataFrame to polrgwas directly.polrgwas(@formula(trait ~ 0 + sex), df, \"../data/hapmap3\")note: Note\nUsers should always make sure that individuals in covariate file or DataFrame match those in Plink fam file. For example, following code checks that the first 2 columns of the covariate.txt file match the first 2 columns of the hapmap3.fam file exactly.covdf = CSV.read(\"../data/covariate.txt\")\nplkfam = CSV.read(\"../data/hapmap3.fam\", header=0, delim=\' \')\nall(covdf[1] .== plkfam[1]) && all(covdf[2] .== plkfam[2])true"
},

{
    "location": "#Timing-1",
    "page": "PolrGWAS.jl",
    "title": "Timing",
    "category": "section",
    "text": "For this moderate-sized data set, polrgwas takes less than 0.2 second.@btime(polrgwas(@formula(trait ~ 0 + sex), \"../data/covariate.txt\", \"../data/hapmap3\"))  128.940 ms (639486 allocations: 29.14 MiB)rm(\"polrgwas.scoretest.txt\")\nrm(\"polrgwas.nullmodel.txt\")"
},

{
    "location": "#Link-functions-1",
    "page": "PolrGWAS.jl",
    "title": "Link functions",
    "category": "section",
    "text": "The link keyword argument of polrgwas can take value LogitLink() (default), ProbitLink() (ordred Probit model), CloglogLink() (proportional hazards model), or CauchyLink().For example, to perform GWAS using the ordred Probit modelpolrgwas(@formula(trait ~ 0 + sex), \"../data/covariate.txt\", \"../data/hapmap3\", \n    link=ProbitLink(), outfile=\"opm\")The estimates in null model and p-values are slightly different from proportional odds moodel.;cat opm.nullmodel.txtStatsModels.DataFrameRegressionModel{PolrModel{Int64,Float64,ProbitLink},Array{Float64,2}}\n\nFormula: trait ~ +sex\n\nCoefficients:\n      Estimate Std.Error  t value Pr(>|t|)\nθ1   -0.866156  0.210746 -4.10995    <1e-4\nθ2   -0.359878   0.20552 -1.75106   0.0809\nθ3    0.247054  0.205135  1.20435   0.2293\nβ1    0.251058  0.128212  1.95814   0.0511;head opm.scoretest.txtchr,pos,snpid,maf,pval\n1,554484,rs10458597,0.0,1.0\n1,758311,rs12562034,0.07763975155279501,0.0064504259191418\n1,967643,rs2710875,0.32407407407407407,1.578504254844719e-5\n1,1168108,rs11260566,0.19158878504672894,4.979075119252628e-6\n1,1375074,rs1312568,0.441358024691358,0.004566574021452637\n1,1588771,rs35154105,0.0,1.0\n1,1789051,rs16824508,0.00462962962962965,0.4819721449163142\n1,1990452,rs2678939,0.4537037037037037,0.33240414006705726\n1,2194615,rs7553178,0.22685185185185186,0.25299085946476113rm(\"opm.nullmodel.txt\")\nrm(\"opm.scoretest.txt\")"
},

{
    "location": "#SNP-and/or-sample-masks-1",
    "page": "PolrGWAS.jl",
    "title": "SNP and/or sample masks",
    "category": "section",
    "text": "In practice, we often perform GWAS on selected SNPs and/or selected samples. They can be specified by the colinds and rowinds keywords of polrgwas function.For example, to perform GWAS on SNPs with minor allele frequency (MAF) above 0.05# create SNP mask\nsnpinds = maf(SnpArray(\"../data/hapmap3.bed\")) .≥ 0.0513928-element BitArray{1}:\n false\n  true\n  true\n  true\n  true\n false\n false\n  true\n  true\n  true\n false\n  true\n  true\n     ⋮\n  true\n  true\n  true\n  true\n  true\n  true\n  true\n  true\n  true\n false\n false\n false# GWAS on selected SNPs\n@time polrgwas(@formula(trait ~ 0 + sex), \"../data/covariate.txt\", \"../data/hapmap3\", \n    colinds = snpinds, outfile=\"commonvariant\")  0.218401 seconds (829.19 k allocations: 39.025 MiB, 3.64% gc time);cat commonvariant.nullmodel.txtStatsModels.DataFrameRegressionModel{PolrModel{Int64,Float64,LogitLink},Array{Float64,2}}\n\nFormula: trait ~ +sex\n\nCoefficients:\n      Estimate Std.Error  t value Pr(>|t|)\nθ1    -1.48564  0.358713 -4.14157    <1e-4\nθ2   -0.569479  0.340649 -1.67175   0.0956\nθ3    0.429815  0.339266   1.2669   0.2061\nβ1    0.424656  0.213911   1.9852   0.0480;head -20 commonvariant.scoretest.txtchr,pos,snpid,maf,pval\n1,758311,rs12562034,0.07763975155279501,0.003001230754791795\n1,967643,rs2710875,0.32407407407407407,2.5117214960313727e-5\n1,1168108,rs11260566,0.19158878504672894,1.13731122530895e-5\n1,1375074,rs1312568,0.441358024691358,0.008317358366815325\n1,1990452,rs2678939,0.4537037037037037,0.29988429741739986\n1,2194615,rs7553178,0.22685185185185186,0.16436415589171793\n1,2396747,rs13376356,0.1448598130841121,0.5372089713885622\n1,2823603,rs1563468,0.4830246913580247,0.23123684490822347\n1,3025087,rs6690373,0.2538699690402477,0.700092366400877\n1,3431124,rs12093117,0.1099071207430341,0.4271132018338309\n1,3633945,rs10910017,0.22187500000000004,0.9141485352635613\n1,4096895,rs6702633,0.4752321981424149,0.006373000780228055\n1,4297388,rs684965,0.3055555555555556,0.09402646589417148\n1,4498133,rs11809295,0.0993788819875776,0.0856953578572345\n1,4698713,rs578528,0.32407407407407407,0.06883563182592009\n1,4899946,rs4654471,0.3580246913580247,0.2267196199789606\n1,5100369,rs6681148,0.13157894736842102,0.16154955342712946\n1,5302730,rs10799197,0.4287925696594427,0.6769491555716189\n1,5502779,rs10796400,0.2314814814814815,0.2449957392353961# extra headline in commonvariant.scoretest.txt\ncountlines(\"commonvariant.scoretest.txt\"), count(snpinds)(12086, 12085)rm(\"commonvariant.scoretest.txt\")\nrm(\"commonvariant.nullmodel.txt\")User should be particularly careful when using the rowinds keyword. Selected rows in SnpArray should exactly match the samples in the null model. Otherwise the results are meaningless."
},

{
    "location": "#Likelihood-ratio-test-(LRT)-1",
    "page": "PolrGWAS.jl",
    "title": "Likelihood ratio test (LRT)",
    "category": "section",
    "text": "By default, polrgwas calculates p-value for each SNP using score test. Score test is fast because it doesn\'t require fitting alternative model for each SNP. User can request likelihood ratio test (LRT) using keyword test=:LRT. LRT is much slower but may be more powerful than score test.@time polrgwas(@formula(trait ~ 0 + sex), \"../data/covariate.txt\", \"../data/hapmap3\", \n    test=:LRT, outfile=\"lrt\") 19.977851 seconds (8.79 M allocations: 2.064 GiB, 2.05% gc time)Test result is output to outfile.lrttest.txt file;head -20 lrt.lrttest.txtchr,pos,snpid,maf,effect,pval\n1,554484,rs10458597,0.0,0.0,1.0\n1,758311,rs12562034,0.07763975155279501,-1.0057833719544362,0.0019185836579804134\n1,967643,rs2710875,0.32407407407407407,-0.6488560566295053,1.805050556976133e-5\n1,1168108,rs11260566,0.19158878504672894,-0.9157225669357916,5.8733847126852225e-6\n1,1375074,rs1312568,0.441358024691358,-0.3318136652577252,0.008081022577831297\n1,1588771,rs35154105,0.0,0.0,1.0\n1,1789051,rs16824508,0.00462962962962965,-0.7338026388700327,0.5169027130128576\n1,1990452,rs2678939,0.4537037037037037,-0.1358649923181972,0.29946402200910055\n1,2194615,rs7553178,0.22685185185185186,-0.2512075640440133,0.16151069094439868\n1,2396747,rs13376356,0.1448598130841121,0.12946142026273552,0.5387338201469207\n1,2623158,rs28753913,0.0,0.0,1.0\n1,2823603,rs1563468,0.4830246913580247,0.15515329587697385,0.2312300208157546\n1,3025087,rs6690373,0.2538699690402477,-0.05966638389967735,0.6995722170701131\n1,3225416,rs12043519,0.029320987654321007,1.1761887120778431,0.002016744167744886\n1,3431124,rs12093117,0.1099071207430341,0.18242332995458113,0.43052012312944177\n1,3633945,rs10910017,0.22187500000000004,-0.017935692627049582,0.9142024828415289\n1,3895935,rs34770924,0.024691358024691357,0.0009448575482692606,0.9980482910119114\n1,4096895,rs6702633,0.4752321981424149,0.4230874102563219,0.0063052493446845445\n1,4297388,rs684965,0.3055555555555556,-0.27091872232225395,0.09226258176075279Note the extra effect column, which is the effect size (regression coefficient) for each SNP. rm(\"lrt.lrttest.txt\")\nrm(\"lrt.nullmodel.txt\")In this example, GWAS by score test takes less than 0.2 second, while GWAS by LRT takes about 20 seconds. About 100 fold difference in run time. "
},

{
    "location": "#Score-test-for-screening,-LRT-for-power-1",
    "page": "PolrGWAS.jl",
    "title": "Score test for screening, LRT for power",
    "category": "section",
    "text": "For large data sets, a practical solution is to perform score test first, then re-do LRT for the most promising SNPs according to score test p-values.Step 1: Perform score test GWAS, results in hapmap.scoretest.txt.@time polrgwas(@formula(trait ~ 0 + sex), \"../data/covariate.txt\", \"../data/hapmap3\", \n    test=:score, outfile=\"hapmap\", verbose=false)  0.233052 seconds (717.09 k allocations: 33.225 MiB, 10.99% gc time);head -20 hapmap.scoretest.txtchr,pos,snpid,maf,pval\n1,554484,rs10458597,0.0,1.0\n1,758311,rs12562034,0.07763975155279501,0.003001230754791795\n1,967643,rs2710875,0.32407407407407407,2.5117214960313727e-5\n1,1168108,rs11260566,0.19158878504672894,1.13731122530895e-5\n1,1375074,rs1312568,0.441358024691358,0.008317358366815325\n1,1588771,rs35154105,0.0,1.0\n1,1789051,rs16824508,0.00462962962962965,0.5274428530031774\n1,1990452,rs2678939,0.4537037037037037,0.29988429741739986\n1,2194615,rs7553178,0.22685185185185186,0.16436415589171793\n1,2396747,rs13376356,0.1448598130841121,0.5372089713885622\n1,2623158,rs28753913,0.0,1.0\n1,2823603,rs1563468,0.4830246913580247,0.23123684490822347\n1,3025087,rs6690373,0.2538699690402477,0.700092366400877\n1,3225416,rs12043519,0.029320987654321007,0.001067589332379961\n1,3431124,rs12093117,0.1099071207430341,0.4271132018338309\n1,3633945,rs10910017,0.22187500000000004,0.9141485352635613\n1,3895935,rs34770924,0.024691358024691357,0.9989137356071247\n1,4096895,rs6702633,0.4752321981424149,0.006373000780228055\n1,4297388,rs684965,0.3055555555555556,0.09402646589417148Step 2: Sort score test p-values and find top 10 SNPs.scorepvals = CSV.read(\"hapmap.scoretest.txt\")[5]13928-element Array{Union{Missing, Float64},1}:\n 1.0                  \n 0.003001230754791795 \n 2.5117214960313727e-5\n 1.13731122530895e-5  \n 0.008317358366815325 \n 1.0                  \n 0.5274428530031774   \n 0.29988429741739986  \n 0.16436415589171793  \n 0.5372089713885622   \n 1.0                  \n 0.23123684490822347  \n 0.700092366400877    \n ⋮                    \n 0.4533525658437616   \n 0.6409263039143044   \n 0.15819888818088343  \n 0.5643246521154451   \n 0.5016172691333742   \n 0.16467049275061096  \n 0.6977871465129508   \n 0.9120338117720298   \n 0.6977871465129944   \n 0.3136740371197593   \n 0.2477983233224063   \n 0.24848633880732612tophits = sortperm(scorepvals)[1:10]10-element Array{Int64,1}:\n  6063\n 13544\n  5071\n  2458\n     4\n  3291\n 13737\n  4183\n     3\n 12727scorepvals[tophits]10-element Array{Union{Missing, Float64},1}:\n 1.4416637280589482e-6\n 5.127005866008025e-6 \n 5.910845635149386e-6 \n 8.756466364522906e-6 \n 1.13731122530895e-5  \n 1.7228784857112124e-5\n 2.008908938410486e-5 \n 2.283909776161149e-5 \n 2.5117214960313727e-5\n 2.949378996752777e-5Step 3: Re-do LRT on top hits.@time polrgwas(@formula(trait ~ 0 + sex), \"../data/covariate.txt\", \"../data/hapmap3\", \n    colinds=tophits, test=:LRT, outfile=\"hapmap\")  0.158214 seconds (408.84 k allocations: 21.313 MiB, 4.23% gc time);head -20 hapmap.lrttest.txtchr,pos,snpid,maf,effect,pval\n1,967643,rs2710875,0.32407407407407407,-0.6488560566295053,1.805050556976133e-5\n1,1168108,rs11260566,0.19158878504672894,-0.9157225669357916,5.8733847126852225e-6\n3,36821790,rs4678553,0.23456790123456794,0.7424952268973525,1.1303825016262592e-5\n4,11017683,rs16881446,0.27554179566563464,-0.7870581482955515,1.1105427468799613e-5\n5,3739190,rs12521166,0.0679012345679012,1.1468852997925303,4.7812882296576845e-5\n6,7574576,rs1885466,0.17746913580246915,0.8750621092263018,7.272346896740631e-6\n7,41152376,rs28880,0.3379629629629629,-0.8146339024453509,9.180126530294943e-7\n20,39225952,rs2076145,0.04475308641975306,1.4412437719976467,5.94831595593157e-5\n23,81247423,rs5923282,0.0030864197530864335,24.22432454878159,0.0016846467294181053\n23,121048059,rs1937165,0.4380804953560371,0.5392313636256603,1.9754643855527356e-5rm(\"hapmap.nullmodel.txt\")\nrm(\"hapmap.lrttest.txt\")\nrm(\"hapmap.scoretest.txt\")"
},

{
    "location": "#GxE-or-other-interactions-1",
    "page": "PolrGWAS.jl",
    "title": "GxE or other interactions",
    "category": "section",
    "text": "In many applications, we want to test SNP effect and/or its interaction with other terms. testformula keyword specifies the test unit besides the covariates in nullformula. In following example, keyword testformula=@formula(trait ~ 0 + snp + snp & sex) instructs polrgwas to test joint effect of snp and snp & sex interaction.polrgwas(@formula(trait ~ 0 + sex), \"../data/covariate.txt\", \"../data/hapmap3\", \n    outfile=\"GxE\", testformula=@formula(trait ~ 0 + snp + snp & sex));head -20 GxE.scoretest.txtchr,pos,snpid,maf,pval\n1,554484,rs10458597,0.0,1.0\n1,758311,rs12562034,0.07763975155279501,0.011788370857064464\n1,967643,rs2710875,0.32407407407407407,0.0001389561913412981\n1,1168108,rs11260566,0.19158878504672894,4.4462425016377604e-5\n1,1375074,rs1312568,0.441358024691358,0.029488132169322546\n1,1588771,rs35154105,0.0,1.0\n1,1789051,rs16824508,0.00462962962962965,0.2922197942259074\n1,1990452,rs2678939,0.4537037037037037,0.3766426135602989\n1,2194615,rs7553178,0.22685185185185186,0.30692923752492324\n1,2396747,rs13376356,0.1448598130841121,0.8215334624145896\n1,2623158,rs28753913,0.0,1.0\n1,2823603,rs1563468,0.4830246913580247,0.3959379791900657\n1,3025087,rs6690373,0.2538699690402477,0.9284920594746979\n1,3225416,rs12043519,0.029320987654321007,0.0046530294043420975\n1,3431124,rs12093117,0.1099071207430341,0.47569653308676063\n1,3633945,rs10910017,0.22187500000000004,0.9880330237520744\n1,3895935,rs34770924,0.024691358024691357,0.3750554140393971\n1,4096895,rs6702633,0.4752321981424149,0.012651490097329717\n1,4297388,rs684965,0.3055555555555556,0.1477240111436921rm(\"GxE.nullmodel.txt\")\nrm(\"GxE.scoretest.txt\")To test the linear and quadratic SNP effects jointlypolrgwas(@formula(trait ~ 0 + sex), \"../data/covariate.txt\", \"../data/hapmap3\", \n    outfile=\"quadratic\", testformula=@formula(trait ~ 0 + snp + snp & snp));head -20 quadratic.scoretest.txtchr,pos,snpid,maf,pval\n1,554484,rs10458597,0.0,1.0\n1,758311,rs12562034,0.07763975155279501,0.012236157097696293\n1,967643,rs2710875,0.32407407407407407,1.0\n1,1168108,rs11260566,0.19158878504672894,6.553991614163971e-5\n1,1375074,rs1312568,0.441358024691358,1.0\n1,1588771,rs35154105,0.0,1.0\n1,1789051,rs16824508,0.00462962962962965,1.0\n1,1990452,rs2678939,0.4537037037037037,1.0\n1,2194615,rs7553178,0.22685185185185186,1.0\n1,2396747,rs13376356,0.1448598130841121,1.0\n1,2623158,rs28753913,0.0,1.0\n1,2823603,rs1563468,0.4830246913580247,1.0\n1,3025087,rs6690373,0.2538699690402477,0.9284972163883947\n1,3225416,rs12043519,0.029320987654321007,1.0\n1,3431124,rs12093117,0.1099071207430341,1.0\n1,3633945,rs10910017,0.22187500000000004,1.0\n1,3895935,rs34770924,0.024691358024691357,1.0\n1,4096895,rs6702633,0.4752321981424149,1.0\n1,4297388,rs684965,0.3055555555555556,1.0rm(\"quadratic.nullmodel.txt\")\nrm(\"quadratic.scoretest.txt\")"
},

{
    "location": "#Docker-1",
    "page": "PolrGWAS.jl",
    "title": "Docker",
    "category": "section",
    "text": "For ease of using PolrGWAS, we provide a Dockerfile so users don\'t need to install Julia and required packages. Only Docker app needs to be installed. Following is tested on ???"
},

{
    "location": "#Plots-1",
    "page": "PolrGWAS.jl",
    "title": "Plots",
    "category": "section",
    "text": "To plot the GWAS results, use the MendelPlots(???) package."
},

]}
