using PolrGWAS, Test, CSV, SnpArrays

const datadir = joinpath(dirname(@__FILE__), "..", "data")
const covfile = datadir * "/covariate.txt"
const plkfile = datadir * "/hapmap3"

@testset "score test" begin
    @time polrgwas(@formula(trait ~ 0 + sex), covfile, plkfile, test=:score)
    @test isfile("polrgwas.nullmodel.txt")
    @test isfile("polrgwas.scoretest.txt")
    scorepvals = CSV.read("polrgwas.scoretest.txt")[5][1:5]
    @test isapprox(scorepvals, [1.0, 4.56531284e-3, 3.10828383e-5, 1.21686724e-5, 8.20686005e-3], rtol=1e-4)
    rm("polrgwas.nullmodel.txt")
    rm("polrgwas.scoretest.txt")
end

@testset "LRT test" begin
    @time polrgwas(@formula(trait ~ 0 + sex), covfile, plkfile, test=:LRT)
    @test isfile("polrgwas.nullmodel.txt")
    @test isfile("polrgwas.lrttest.txt")
    lrtpvals = CSV.read("polrgwas.lrttest.txt")[6][1:5]
    @test isapprox(lrtpvals, [1.0, 1.91858366e-3, 1.80505056e-5, 5.87338471e-6, 8.08102258e-3], rtol=1e-4)
    rm("polrgwas.nullmodel.txt")
    rm("polrgwas.lrttest.txt")
end

@testset "snpmodel" begin
    # dominant model 
    polrgwas(@formula(trait ~ 0 + sex), covfile, plkfile, test=:score, snpmodel=DOMINANT_MODEL)
    @test isfile("polrgwas.nullmodel.txt")
    @test isfile("polrgwas.scoretest.txt")
    scorepvals = CSV.read("polrgwas.scoretest.txt")[5][1:5]
    @test isapprox(scorepvals, [1.0, 0.14295, 0.000471942, 0.00555348, 0.000652844], rtol=1e-4)
    rm("polrgwas.nullmodel.txt")
    rm("polrgwas.scoretest.txt")
    # recessive model 
    polrgwas(@formula(trait ~ 0 + sex), covfile, plkfile, test=:score, snpmodel=RECESSIVE_MODEL)
    @test isfile("polrgwas.nullmodel.txt")
    @test isfile("polrgwas.scoretest.txt")
    scorepvals = CSV.read("polrgwas.scoretest.txt")[5][1:5]
    @test isapprox(scorepvals, [1.0, 0.00673612, 0.000279908, 4.15322e-5, 0.167642], rtol=1e-4)
    rm("polrgwas.nullmodel.txt")
    rm("polrgwas.scoretest.txt")
end

@testset "link" begin
    polrgwas(@formula(trait ~ 0 + sex), covfile, plkfile, link=ProbitLink(), outfile="opm")
    @test isfile("opm.nullmodel.txt")
    @test isfile("opm.scoretest.txt")
    scorepvals = CSV.read("opm.scoretest.txt")[5][1:5]
    @test isapprox(scorepvals, [1.0, 1.00769167e-2, 2.62725649e-5, 1.08974849e-5, 5.10288399e-3], rtol=1e-4)
    rm("opm.nullmodel.txt")
    rm("opm.scoretest.txt")
end

@testset "mask" begin
    polrgwas(@formula(trait ~ 0 + sex), covfile, plkfile, colinds=1:5, outfile="first5snps")
    @test isfile("first5snps.nullmodel.txt")
    @test isfile("first5snps.scoretest.txt")
    @test countlines("first5snps.scoretest.txt") == 6
    scorepvals = CSV.read("first5snps.scoretest.txt")[5]
    @test isapprox(scorepvals, [1.0, 4.56531284e-3, 3.10828383e-5, 1.21686724e-5, 8.20686005e-3], rtol=1e-4)
    rm("first5snps.nullmodel.txt")
    rm("first5snps.scoretest.txt")
end

@testset "test formula" begin
    # score test
    polrgwas(@formula(trait ~ 0 + sex), covfile, plkfile, outfile="GxE", testformula=@formula(trait ~ 0 + snp + snp & sex))
    @test isfile("GxE.nullmodel.txt")
    @test isfile("GxE.scoretest.txt")
    scorepvals = CSV.read("GxE.scoretest.txt")[5][1:5]
    @test isapprox(scorepvals, [1.0, 1.74460104e-2, 1.66707324e-4, 4.76376246e-5, 2.91384712e-2], rtol=1e-4)
    rm("GxE.nullmodel.txt")
    rm("GxE.scoretest.txt")
    # LRT, only first 5 SNPs
    polrgwas(@formula(trait ~ 0 + sex), covfile, plkfile, outfile="GxE", 
    testformula=@formula(trait ~ 0 + snp + snp & sex), test=:LRT, colinds=1:5)
    @test isfile("GxE.nullmodel.txt")
    @test isfile("GxE.lrttest.txt")
    lrtpvals = CSV.read("GxE.lrttest.txt")[end]
    @show lrtpvals
    @test isapprox(lrtpvals, [1.0, 7.22410973e-3, 1.01730983e-4, 1.88174211e-5, 2.88295705e-2], rtol=1e-4)
    rm("GxE.nullmodel.txt")
    rm("GxE.lrttest.txt")
end
