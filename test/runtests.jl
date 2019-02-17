using OrdinalGWAS, Test, CSV, SnpArrays

const datadir = joinpath(dirname(@__FILE__), "..", "data")
const covfile = datadir * "/covariate.txt"
const plkfile = datadir * "/hapmap3"

@testset "score test" begin
    @time ordinalgwas(@formula(trait ~ sex), covfile, plkfile, test=:score)
    @test isfile("ordinalgwas.null.txt")
    @test isfile("ordinalgwas.pval.txt")
    scorepvals = open(CSV.read, "ordinalgwas.pval.txt")[5][1:5]
    @test isapprox(scorepvals, [1.0, 4.56531284e-3, 3.10828383e-5, 1.21686724e-5, 8.20686005e-3], rtol=1e-4)
    rm("ordinalgwas.null.txt")
    rm("ordinalgwas.pval.txt")
end

@testset "LRT test" begin
    @time ordinalgwas(@formula(trait ~ sex), covfile, plkfile, test=:LRT)
    @test isfile("ordinalgwas.null.txt")
    @test isfile("ordinalgwas.pval.txt")
    lrtpvals = open(CSV.read, "ordinalgwas.pval.txt")[6][1:5]
    @test isapprox(lrtpvals, [1.0, 1.91858366e-3, 1.80505056e-5, 5.87338471e-6, 8.08102258e-3], rtol=1e-4)
    rm("ordinalgwas.null.txt")
    rm("ordinalgwas.pval.txt")
end

@testset "snpmodel" begin
    # dominant model 
    ordinalgwas(@formula(trait ~ sex), covfile, plkfile, test=:score, snpmodel=DOMINANT_MODEL)
    @test isfile("ordinalgwas.null.txt")
    @test isfile("ordinalgwas.pval.txt")
    scorepvals = open(CSV.read, "ordinalgwas.pval.txt")[5][1:5]
    @test isapprox(scorepvals, [1.0, 0.14295, 0.000471942, 0.00555348, 0.000652844], rtol=1e-4)
    rm("ordinalgwas.null.txt")
    rm("ordinalgwas.pval.txt")
    # recessive model 
    ordinalgwas(@formula(trait ~ sex), covfile, plkfile, test=:score, snpmodel=RECESSIVE_MODEL)
    @test isfile("ordinalgwas.null.txt")
    @test isfile("ordinalgwas.pval.txt")
    scorepvals = open(CSV.read, "ordinalgwas.pval.txt")[5][1:5]
    @test isapprox(scorepvals, [1.0, 0.00673612, 0.000279908, 4.15322e-5, 0.167642], rtol=1e-4)
    rm("ordinalgwas.null.txt")
    rm("ordinalgwas.pval.txt")
end

@testset "link" begin
    ordinalgwas(@formula(trait ~ sex), covfile, plkfile, link=ProbitLink(), pvalfile="opm.pval.txt")
    @test isfile("ordinalgwas.null.txt")
    @test isfile("opm.pval.txt")
    scorepvals = open(CSV.read, "opm.pval.txt")[5][1:5]
    @test isapprox(scorepvals, [1.0, 1.00769167e-2, 2.62725649e-5, 1.08974849e-5, 5.10288399e-3], rtol=1e-4)
    rm("ordinalgwas.null.txt")
    rm("opm.pval.txt")
end

@testset "snp mask" begin
    ordinalgwas(@formula(trait ~ sex), covfile, plkfile, snpinds=1:5, pvalfile="first5snps.pval.txt")
    @test isfile("ordinalgwas.null.txt")
    @test isfile("first5snps.pval.txt")
    @test countlines("first5snps.pval.txt") == 6
    scorepvals = open(CSV.read, "first5snps.pval.txt")[5]
    @test isapprox(scorepvals, [1.0, 4.56531284e-3, 3.10828383e-5, 1.21686724e-5, 8.20686005e-3], rtol=1e-4)
    rm("ordinalgwas.null.txt")
    rm("first5snps.pval.txt")
end

@testset "sub samples" begin
    # only use first 300 samples
    @time ordinalgwas(@formula(trait ~ sex), covfile, plkfile, test=:score, covrowinds=1:300, bedrowinds=1:300)
    @test isfile("ordinalgwas.null.txt")
    @test isfile("ordinalgwas.pval.txt")
    scorepvals = open(CSV.read, "ordinalgwas.pval.txt")[5][1:5]
    @test isapprox(scorepvals, [1.0, 0.00355969, 0.000123604, 5.2213e-6, 0.00758234], rtol=1e-4)
    rm("ordinalgwas.null.txt")
    rm("ordinalgwas.pval.txt")
end

@testset "test formula" begin
    # score test
    ordinalgwas(@formula(trait ~ sex), covfile, plkfile, pvalfile="GxE.pval.txt", 
    testformula=@formula(trait ~ snp + snp & sex))
    @test isfile("ordinalgwas.null.txt")
    @test isfile("GxE.pval.txt")
    scorepvals = open(CSV.read, "GxE.pval.txt")[5][1:5]
    @test isapprox(scorepvals, [1.0, 1.74460104e-2, 1.66707324e-4, 4.76376246e-5, 2.91384712e-2], rtol=1e-4)
    rm("ordinalgwas.null.txt")
    rm("GxE.pval.txt")
    # LRT, only first 5 SNPs
    ordinalgwas(@formula(trait ~ sex), covfile, plkfile, pvalfile="GxE.pval.txt", 
    testformula=@formula(trait ~ snp + snp & sex), test=:LRT, snpinds=1:5)
    @test isfile("ordinalgwas.null.txt")
    @test isfile("GxE.pval.txt")
    lrtpvals = open(CSV.read, "GxE.pval.txt")[end]
    @test isapprox(lrtpvals, [1.0, 7.22410973e-3, 1.01730983e-4, 1.88174211e-5, 2.88295705e-2], rtol=1e-4)
    rm("ordinalgwas.null.txt")
    rm("GxE.pval.txt")
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
            pvals_chr1 = open(CSV.read, pvalfile)[5][1:5]
            @test isapprox(pvals_chr1, [1.0, 4.56531284e-3, 3.10828383e-5, 1.21686724e-5, 8.20686005e-3], rtol=1e-4)    
        end
        rm(plinkfile * ".pval.txt")
    end
    # gwas by chromosome, use fitted null model each time, use uncompressed Plink set
    @time for chr in 1:23
        plinkfile = plkfile * ".chr." * string(chr)
        pvalfile = plkfile * ".chr." * string(chr) * ".pval.txt"
        ordinalgwas(nm, plinkfile, pvalfile = pvalfile)
        @test isfile(pvalfile)
        if chr == 1
            pvals_chr1 = open(CSV.read, pvalfile)[5][1:5]
            @test isapprox(pvals_chr1, [1.0, 4.56531284e-3, 3.10828383e-5, 1.21686724e-5, 8.20686005e-3], rtol=1e-4)    
        end
        rm(pvalfile)
    end
    # gwas by chromosome, use fitted null model each time, use compressed bed and bim files
    @time for chr in 1:23
        bedfile = plkfile * ".chr." * string(chr) * ".bed.gz"
        bimfile = plkfile * ".chr." * string(chr) * ".bim.gz"
        pvalfile = plkfile * ".chr." * string(chr) * ".pval.txt"
        ordinalgwas(nm, bedfile, bimfile, 324; pvalfile = pvalfile)
        @test isfile(pvalfile)
        if chr == 1
            pvals_chr1 = open(CSV.read, pvalfile)[5][1:5]
            @test isapprox(pvals_chr1, [1.0, 4.56531284e-3, 3.10828383e-5, 1.21686724e-5, 8.20686005e-3], rtol=1e-4)    
        end
        rm(pvalfile)
    end
    # clean up
    # delete result files
    isfile("ordinalgwas.null.txt") && rm("ordinalgwas.null.txt")
    for chr in 1:26
        plinkfile = plkfile * ".chr." * string(chr)
        # delete uncompressed chromosome Plink files
        isfile(plinkfile * ".bed") && rm(plinkfile * ".bed")
        isfile(plinkfile * ".fam") && rm(plinkfile * ".fam")
        isfile(plinkfile * ".bim") && rm(plinkfile * ".bim")
        # delete compressed chromosome Plink files
        isfile(plinkfile * ".bed.gz") && rm(plinkfile * ".bed.gz")
        isfile(plinkfile * ".fam.gz") && rm(plinkfile * ".fam.gz")
        isfile(plinkfile * ".bim.gz") && rm(plinkfile * ".bim.gz")
        # delete pval files
        isfile(plinkfile * ".pval.txt") && rm(plinkfile * ".pval.txt")
    end
end
