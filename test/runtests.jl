using OrdinalGWAS, Test, CSV, SnpArrays

const datadir = joinpath(dirname(@__FILE__), "..", "data")
const covfile = datadir * "/covariate.txt"
const plkfile = datadir * "/hapmap3"
const snpsetfile = datadir * "/hapmap_snpsetfile.txt"
const vcfcovfile = datadir * "/vcf_example.csv"
const vcffile = datadir * "/vcf_test"
const vcfsnpsetfile = datadir * "/snpsetfile_vcf.txt"


@testset "score test" begin
    @time ordinalgwas(@formula(trait ~ sex), covfile, plkfile, test=:score)
    @test isfile("ordinalgwas.null.txt")
    @test isfile("ordinalgwas.pval.txt")
    scorepvals = open(CSV.read, "ordinalgwas.pval.txt")[!, 6][1:5]
    @test isapprox(scorepvals, [1.0, 4.56531284e-3, 3.10828383e-5, 1.21686724e-5, 8.20686005e-3], rtol=1e-4)
    rm("ordinalgwas.null.txt", force=true)
    rm("ordinalgwas.pval.txt", force=true)

    # VCF
    ordinalgwas(@formula(y ~ sex), vcfcovfile, vcffile; geneticformat = "VCF", 
        vcftype = :DS, geneticrowinds = 1:190, snpinds = [86; 656], 
    test = :score, covrowinds = 1:190)
    scorepvals = open(CSV.read, "ordinalgwas.pval.txt")[!, end][1:2]
    @test isapprox(scorepvals, [0.00762272, 0.000668338], rtol=1e-4)
    rm("ordinalgwas.null.txt", force=true)
    rm("ordinalgwas.pval.txt", force=true)
end

@testset "LRT test" begin
    @time ordinalgwas(@formula(trait ~ sex), covfile, plkfile, test=:LRT)
    @test isfile("ordinalgwas.null.txt")
    @test isfile("ordinalgwas.pval.txt")
    lrtpvals = open(CSV.read, "ordinalgwas.pval.txt")[!, 7][1:5]
    @test isapprox(lrtpvals, [1.0, 1.91858366e-3, 1.80505056e-5, 5.87338471e-6, 8.08102258e-3], rtol=1e-4)
    rm("ordinalgwas.null.txt", force=true)
    rm("ordinalgwas.pval.txt", force=true)

    # VCF
    ordinalgwas(@formula(y ~ sex), vcfcovfile, vcffile; geneticformat = "VCF", 
        vcftype = :GT, snpinds = [86; 656], test = :LRT)
    lrtpvals = open(CSV.read, "ordinalgwas.pval.txt")[!, end][1:2]
    @test isapprox(lrtpvals, [0.00955468405473856, 0.0007086063489553798], rtol=1e-4)
    rm("ordinalgwas.null.txt", force=true)
    rm("ordinalgwas.pval.txt", force=true)
end

@testset "snpmodel" begin
    # dominant model 
    ordinalgwas(@formula(trait ~ sex), covfile, plkfile, test=:score, snpmodel=DOMINANT_MODEL)
    @test isfile("ordinalgwas.null.txt")
    @test isfile("ordinalgwas.pval.txt")
    scorepvals = open(CSV.read, "ordinalgwas.pval.txt")[!, 6][1:5]
    @test isapprox(scorepvals, [1.0, 0.14295, 0.000471942, 0.00555348, 0.000652844], rtol=1e-4)
    rm("ordinalgwas.null.txt", force=true)
    rm("ordinalgwas.pval.txt", force=true)
    # recessive model 
    ordinalgwas(@formula(trait ~ sex), covfile, plkfile, test=:score, snpmodel=RECESSIVE_MODEL)
    @test isfile("ordinalgwas.null.txt")
    @test isfile("ordinalgwas.pval.txt")
    scorepvals = open(CSV.read, "ordinalgwas.pval.txt")[!, 6][1:5]
    @test isapprox(scorepvals, [1.0, 0.00673612, 0.000279908, 4.15322e-5, 0.167642], rtol=1e-4)
    rm("ordinalgwas.null.txt", force=true)
    rm("ordinalgwas.pval.txt", force=true)
end

@testset "link" begin
    ordinalgwas(@formula(trait ~ sex), covfile, plkfile, link=ProbitLink(), pvalfile="opm.pval.txt")
    @test isfile("ordinalgwas.null.txt")
    @test isfile("opm.pval.txt")
    scorepvals = open(CSV.read, "opm.pval.txt")[!, 6][1:5]
    @test isapprox(scorepvals, [1.0, 1.00769167e-2, 2.62725649e-5, 1.08974849e-5, 5.10288399e-3], rtol=1e-4)
    rm("ordinalgwas.null.txt", force=true)
    rm("opm.pval.txt", force=true)
end

@testset "snp mask" begin
    ordinalgwas(@formula(trait ~ sex), covfile, plkfile, snpinds=1:5, pvalfile="first5snps.pval.txt")
    @test isfile("ordinalgwas.null.txt")
    @test isfile("first5snps.pval.txt")
    @test countlines("first5snps.pval.txt") == 6
    scorepvals = open(CSV.read, "first5snps.pval.txt")[!, 6]
    @test isapprox(scorepvals, [1.0, 4.56531284e-3, 3.10828383e-5, 1.21686724e-5, 8.20686005e-3], rtol=1e-4)
    rm("ordinalgwas.null.txt", force=true)
    rm("first5snps.pval.txt", force=true)
end

@testset "sub samples" begin
    # only use first 300 samples
    @time ordinalgwas(@formula(trait ~ sex), covfile, plkfile, test=:score, covrowinds=1:300, geneticrowinds=1:300)
    @test isfile("ordinalgwas.null.txt")
    @test isfile("ordinalgwas.pval.txt")
    scorepvals = open(CSV.read, "ordinalgwas.pval.txt")[!, 6][1:5]
    @test isapprox(scorepvals, [1.0, 0.00355969, 0.000123604, 5.2213e-6, 0.00758234], rtol=1e-4)
    rm("ordinalgwas.null.txt", force=true)
    rm("ordinalgwas.pval.txt", force=true)
end

@testset "test formula" begin
    # score test
    ordinalgwas(@formula(trait ~ sex), covfile, plkfile, pvalfile="GxE.pval.txt", 
        testformula=@formula(trait ~ snp + snp & sex))
    @test isfile("ordinalgwas.null.txt")
    @test isfile("GxE.pval.txt")
    scorepvals = open(CSV.read, "GxE.pval.txt")[!, 6][1:5]
    @test isapprox(scorepvals, [1.0, 1.74460104e-2, 1.66707324e-4, 4.76376246e-5, 2.91384712e-2], rtol=1e-4)
    rm("ordinalgwas.null.txt", force=true)
    rm("GxE.pval.txt", force=true)
    # LRT, only first 5 SNPs
    ordinalgwas(@formula(trait ~ sex), covfile, plkfile, pvalfile="GxE.pval.txt", 
        testformula=@formula(trait ~ snp + snp & sex), test=:LRT, snpinds=1:5)
    @test isfile("ordinalgwas.null.txt")
    @test isfile("GxE.pval.txt")
    lrtpvals = open(CSV.read, "GxE.pval.txt")[!, end]
    @test isapprox(lrtpvals, [1.0, 7.22410973e-3, 1.01730983e-4, 1.88174211e-5, 2.88295705e-2], rtol=1e-4)
    rm("ordinalgwas.null.txt", force=true)
    rm("GxE.pval.txt", force=true)
end

@testset "snpset" begin
    #window
    #score test
    ordinalgwas(@formula(trait ~ sex), covfile, plkfile, pvalfile = "snpset.pval.txt",
        snpset=250, analysistype = "snpset")
    @test isfile("ordinalgwas.null.txt")
    @test isfile("snpset.pval.txt")
    scorepvals = open(CSV.read, "snpset.pval.txt")[!, end][1:5]
    #@test isapprox(scorepvals, [1.0, 1.74460104e-2, 1.66707324e-4, 4.76376246e-5, 2.91384712e-2], rtol=1e-4)
    rm("ordinalgwas.null.txt", force=true)
    rm("snpset.pval.txt", force=true)
    #lrt 
    ordinalgwas(@formula(trait ~ sex), covfile, plkfile, pvalfile = "snpset.pval.txt",
        snpset=25, test=:LRT, analysistype = "snpset")
    @test isfile("ordinalgwas.null.txt")
    @test isfile("snpset.pval.txt")
    lrtpvals = open(CSV.read, "snpset.pval.txt")[!, end][1:5]
    @test isapprox(lrtpvals, [2.1817554071810948e-13, 0.2865769729670889, 0.32507802233937966,
    0.3344823237332578, 0.42948375949508427], rtol=1e-4)
    rm("ordinalgwas.null.txt", force=true)
    rm("snpset.pval.txt", force=true)

    # VCF
    ordinalgwas(@formula(y ~ sex), vcfcovfile, vcffile; geneticformat = "VCF", 
        vcftype = :DS, pvalfile = "snpset.pval.txt",
    snpset=250, test=:score, analysistype = "snpset")
    @test isfile("ordinalgwas.null.txt")
    @test isfile("snpset.pval.txt")
    scorepvals = open(CSV.read, "snpset.pval.txt")[!, end][1:5]
    @test isapprox(scorepvals, [0.4278366599084349, 
    0.42781616067453476, 0.4519573757701432, 
    0.4278763804444088, 0.4345883185481474], rtol=1e-4)
    rm("ordinalgwas.null.txt", force=true)
    rm("snpset.pval.txt", force=true)

    ordinalgwas(@formula(y ~ sex), vcfcovfile, vcffile; geneticformat = "VCF", 
        vcftype = :DS, pvalfile = "snpset.pval.txt",
    snpset=25, test=:LRT, analysistype = "snpset")
    lrtpvals = open(CSV.read, "snpset.pval.txt")[!, end][1:5]
    @test isapprox(lrtpvals, [1.0
    0.9999996378252629
    1.0
    0.9999999999996334
    0.9999999976994737], rtol=1e-4)
    rm("ordinalgwas.null.txt", force=true)
    rm("snpset.pval.txt", force=true)

    #snpset file
    #score test
    ordinalgwas(@formula(trait ~ sex), covfile, plkfile, pvalfile = "snpset.pval.txt",
        snpset = snpsetfile, analysistype = "snpset")
    @test isfile("ordinalgwas.null.txt")
    @test isfile("snpset.pval.txt")
    scorepvals = open(CSV.read, "snpset.pval.txt")[!, end][1:5]
    @test isapprox(scorepvals, [1.72134e-5, 0.036925, 0.747855,
     0.0276508, 0.611958], rtol=1e-4)
    rm("ordinalgwas.null.txt", force=true)
    rm("snpset.pval.txt", force=true)
    #lrt 
    ordinalgwas(@formula(trait ~ sex), covfile, plkfile, pvalfile = "snpset.pval.txt",
        snpset = snpsetfile, test = :lrt, analysistype = "snpset")
    @test isfile("ordinalgwas.null.txt")
    @test isfile("snpset.pval.txt")
    lrtpvals = open(CSV.read, "snpset.pval.txt")[!, end][1:5]
    @test isapprox(lrtpvals, [6.75377e-13, 0.000256566, 0.359382,
     0.000163268, 0.0867508], rtol=1e-4)
    rm("ordinalgwas.null.txt", force=true)
    rm("snpset.pval.txt", force=true)

    # VCF
    ordinalgwas(@formula(y ~ sex), vcfcovfile, vcffile; geneticformat = "VCF", 
        vcftype = :DS, pvalfile = "snpset.pval.txt",
        snpset = vcfsnpsetfile, test = :score, analysistype = "snpset")
    @test isfile("ordinalgwas.null.txt")
    @test isfile("snpset.pval.txt")
    scorepvals = open(CSV.read, "snpset.pval.txt")[!, end][1:5]
    @test isapprox(scorepvals, [0.06814002639277685, 0.5566664123188036, 
    0.520381855174413, 0.07557764137466122, 0.5620803022597403], rtol=1e-4)
    rm("ordinalgwas.null.txt", force=true)
    rm("snpset.pval.txt", force=true)

    ordinalgwas(@formula(y ~ sex), vcfcovfile, vcffile; geneticformat = "VCF", 
        vcftype = :DS, pvalfile = "snpset.pval.txt",
        snpset = vcfsnpsetfile, test = :lrt, analysistype = "snpset")
    @test isfile("ordinalgwas.null.txt")
    @test isfile("snpset.pval.txt")
    lrtpvals = open(CSV.read, "snpset.pval.txt")[!, end][1:5]
    @test isapprox(lrtpvals, [0.09069975735216675, 0.6465153355309161, 
    0.6307411986741357, 0.06275888993714969, 0.50252192003468], rtol=1e-4)
    rm("ordinalgwas.null.txt", force=true)
    rm("snpset.pval.txt", force=true)

    #specific snp (one snpset)
    #score test
    ordinalgwas(@formula(trait ~ sex), covfile, plkfile, pvalfile = "snpset.pval.txt",
        snpset = 50:55, analysistype = "snpset")
    @test isfile("ordinalgwas.null.txt")
    @test isfile("snpset.pval.txt")
    scorepvals = open("snpset.pval.txt")
    scorepval = split(readline(scorepvals))[end]
    close(scorepvals)
    @test isapprox(parse(Float64, scorepval), 0.3647126536663949, rtol=1e-4)
    rm("ordinalgwas.null.txt", force=true)
    rm("snpset.pval.txt", force=true)
    #lrt 
    ordinalgwas(@formula(trait ~ sex), covfile, plkfile, pvalfile = "snpset.pval.txt",
        snpset = collect(1:15), test=:LRT, analysistype = "snpset")
    @test isfile("ordinalgwas.null.txt")
    @test isfile("snpset.pval.txt")
    lrtpvals = open("snpset.pval.txt")
    lrtpval = split(readline(lrtpvals))[end]
    close(lrtpvals)
    @test isapprox(parse(Float64, lrtpval), 7.525696044086955e-15, rtol=1e-4)
    rm("ordinalgwas.null.txt", force=true)
    rm("snpset.pval.txt", force=true)

    # VCF
    ordinalgwas(@formula(y ~ sex), vcfcovfile, vcffile; geneticformat = "VCF", 
        vcftype = :DS, pvalfile = "snpset.pval.txt",
        snpset=85:90, test=:score, analysistype = "snpset")
    @test isfile("ordinalgwas.null.txt")
    @test isfile("snpset.pval.txt")
    scorepvals = open("snpset.pval.txt")
    scorepval = split(readline(scorepvals))[end]
    close(scorepvals)
    @test isapprox(parse(Float64, scorepval), 0.0965927460813927, rtol=1e-4)

    ordinalgwas(@formula(y ~ sex), vcfcovfile, vcffile; geneticformat = "VCF", 
        vcftype = :DS, pvalfile = "snpset.pval.txt",
        snpset=85:90, test=:lrt, analysistype = "snpset")
    lrtpvals = open("snpset.pval.txt")
    lrtpval = split(readline(lrtpvals))[end]
    close(lrtpvals)
    @test isapprox(parse(Float64, lrtpval), 0.0732485446883825, rtol=1e-4)
    rm("ordinalgwas.null.txt", force=true)
    rm("snpset.pval.txt", force=true)
end

@testset "GxE snp in null" begin
    ordinalgwas(@formula(trait ~ sex), covfile, plkfile, e = :sex, pvalfile = "gxe_snp.pval.txt",
        snpinds=1:5, test=:score, analysistype = "gxe")
    @test isfile("gxe_snp.pval.txt")
    scorepvals = open(CSV.read, "gxe_snp.pval.txt")[!, end][1:5]
    @test isapprox(scorepvals, [1.0, 0.637742242597749, 0.9667114198051628,
    0.26352674694121003, 0.7811133315582837], rtol=1e-4)
    rm("gxe_snp.pval.txt", force=true)

    ordinalgwas(@formula(trait ~ sex), covfile, plkfile, e = "sex", pvalfile = "gxe_snp.pval.txt",
        snpinds=1:5, test=:LRT, analysistype = "gxe")
    @test isfile("gxe_snp.pval.txt")
    lrtpvals = open(CSV.read, "gxe_snp.pval.txt")[!, end][1:5]
    @test isapprox(lrtpvals, [1.0, 0.6279730133445315, 0.9671662821946985,
    0.26693502209463904, 0.7810214899265426], rtol=1e-4)
    rm("gxe_snp.pval.txt", force=true)

    # VCF
    ordinalgwas(@formula(y ~ sex), vcfcovfile, vcffile, e = :sex; geneticformat = "VCF", 
        vcftype = :DS, pvalfile = "gxe_snp.pval.txt",
        snpinds=1:5, test=:score, analysistype = "gxe")
    scorepvals = open(CSV.read, "gxe_snp.pval.txt")[!, end][1:5]
    @test isapprox(scorepvals, [0.45861769035708144, 1.0, 0.03804677528312195,
     0.18254151103030725, 0.34454453512541156], rtol=1e-4)
    rm("gxe_snp.pval.txt", force=true)

    ordinalgwas(@formula(y ~ sex), vcfcovfile, vcffile, e = :sex; geneticformat = "VCF", 
        vcftype = :DS, pvalfile = "gxe_snp.pval.txt",
        snpinds=1:5, test=:lrt, analysistype = "gxe")
    @test isfile("gxe_snp.pval.txt")
    lrtpvals = open(CSV.read, "gxe_snp.pval.txt")[!, end][1:5]
    @test isapprox(lrtpvals, [0.526667096902957, 1.0, 0.008073040021982156, 
    0.10590569987122991, 0.3557829099471382], rtol=1e-4)
    rm("gxe_snp.pval.txt", force=true)
end

@testset "split, gz" begin
    # split hapmap3 by chromosome
    SnpArrays.split_plink(plkfile, :chromosome; prefix = datadir * "/hapmap3.chr.")
    # compress to gz
    for chr in 1:23
        plinkfile = plkfile * ".chr." * string(chr)
        SnpArrays.compress_plink(plinkfile)
        @test isfile(plinkfile * ".bed.gz")
        @test isfile(plinkfile * ".fam.gz")
        @test isfile(plinkfile * ".bim.gz")
    end
    # fit null model
    @time nm = ordinalgwas(@formula(trait ~ sex), covfile, nothing)
    @test isfile("ordinalgwas.null.txt")
    # gwas by chromosome, refit null model each time, use uncompressed Plink set
    @time for chr in 1:23
        plinkfile = plkfile * ".chr." * string(chr)
        pvalfile = plkfile * ".chr." * string(chr) * ".pval.txt"
        ordinalgwas(@formula(trait ~ sex), covfile, plinkfile, pvalfile = pvalfile)
        @test isfile(pvalfile)
        if chr == 1
            pvals_chr1 = open(CSV.read, pvalfile)[!, 6][1:5]
            @test isapprox(pvals_chr1, [1.0, 4.56531284e-3, 3.10828383e-5, 1.21686724e-5, 8.20686005e-3], rtol=1e-4)    
        end
        rm(plinkfile * ".pval.txt", force=true)
    end
    # gwas by chromosome, use fitted null model each time, use uncompressed Plink set
    @time for chr in 1:23
        plinkfile = plkfile * ".chr." * string(chr)
        pvalfile = plkfile * ".chr." * string(chr) * ".pval.txt"
        ordinalgwas(nm, plinkfile, pvalfile = pvalfile)
        @test isfile(pvalfile)
        if chr == 1
            pvals_chr1 = open(CSV.read, pvalfile)[!, 6][1:5]
            @test isapprox(pvals_chr1, [1.0, 4.56531284e-3, 3.10828383e-5, 1.21686724e-5, 8.20686005e-3], rtol=1e-4)    
        end
        rm(pvalfile, force=true)
    end
    # gwas by chromosome, use fitted null model each time, use compressed bed and bim files
    @time for chr in 1:23
        bedfile = plkfile * ".chr." * string(chr) * ".bed.gz"
        bimfile = plkfile * ".chr." * string(chr) * ".bim.gz"
        pvalfile = plkfile * ".chr." * string(chr) * ".pval.txt"
        ordinalgwas(nm, bedfile, bimfile, 324; pvalfile = pvalfile)
        @test isfile(pvalfile)
        if chr == 1
            pvals_chr1 = open(CSV.read, pvalfile)[!, 6][1:5]
            @test isapprox(pvals_chr1, [1.0, 4.56531284e-3, 3.10828383e-5, 1.21686724e-5, 8.20686005e-3], rtol=1e-4)    
        end
        rm(pvalfile, force=true)
    end
    # clean up
    # delete result files
    isfile("ordinalgwas.null.txt") && rm("ordinalgwas.null.txt")
    for chr in 1:26
        plinkfile = plkfile * ".chr." * string(chr)
        # delete uncompressed chromosome Plink files
        rm(plinkfile * ".bed", force=true)
        rm(plinkfile * ".fam", force=true)
        rm(plinkfile * ".bim", force=true)
        # delete compressed chromosome Plink files
        rm(plinkfile * ".bed.gz", force=true)
        rm(plinkfile * ".fam.gz", force=true)
        rm(plinkfile * ".bim.gz", force=true)
        # delete pval files
        rm(plinkfile * ".pval.txt", force=true)
    end
end
