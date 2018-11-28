using PolrGWAS, Test, CSV, SnpArrays

const datadir = joinpath(dirname(@__FILE__), "..", "data")
const covfile = datadir * "/covariate.txt"
const plkfile = datadir * "/hapmap3"

@testset "score test" begin
    @time polrgwas(@formula(trait ~ 0 + sex), covfile, plkfile, test=:score)
    @test isfile("polrgwas.nullmodel.txt")
    @test isfile("polrgwas.scoretest.txt")
    scorepvals = CSV.read("polrgwas.scoretest.txt")[5][1:5]
    @test all(scorepvals .≈ [1.0, 3.00123075e-3, 2.51172149e-5, 1.137311225e-5, 8.31735837e-3])
    rm("polrgwas.nullmodel.txt")
    rm("polrgwas.scoretest.txt")
end

@testset "LRT test" begin
    @time polrgwas(@formula(trait ~ 0 + sex), covfile, plkfile, test=:LRT)
    @test isfile("polrgwas.nullmodel.txt")
    @test isfile("polrgwas.lrttest.txt")
    lrtpvals = CSV.read("polrgwas.lrttest.txt")[6][1:5]
    @test all(lrtpvals .≈ [1.0, 1.91858366e-3, 1.80505056e-5, 5.87338471e-6, 8.08102258e-3])
    rm("polrgwas.nullmodel.txt")
    rm("polrgwas.lrttest.txt")
end

@testset "snpmodel" begin
    # dominant model 
    polrgwas(@formula(trait ~ 0 + sex), covfile, plkfile, test=:score, snpmodel=DOMINANT_MODEL)
    @test isfile("polrgwas.nullmodel.txt")
    @test isfile("polrgwas.scoretest.txt")
    scorepvals = CSV.read("polrgwas.scoretest.txt")[5][1:5]
    @test all(scorepvals .≈ [1.0, 1.24849027e-1, 4.08779119e-4, 5.19281523e-3, 7.04007765e-4])
    rm("polrgwas.nullmodel.txt")
    rm("polrgwas.scoretest.txt")
    # recessive model 
    polrgwas(@formula(trait ~ 0 + sex), covfile, plkfile, test=:score, snpmodel=RECESSIVE_MODEL)
    @test isfile("polrgwas.nullmodel.txt")
    @test isfile("polrgwas.scoretest.txt")
    scorepvals = CSV.read("polrgwas.scoretest.txt")[5][1:5]
    @test all(scorepvals .≈ [1.0, 4.56151864e-3, 2.77072395e-4, 4.25357392e-5, 1.66390827e-1])
    rm("polrgwas.nullmodel.txt")
    rm("polrgwas.scoretest.txt")
end

@testset "link" begin
    polrgwas(@formula(trait ~ 0 + sex), covfile, plkfile, link=ProbitLink(), outfile="opm")
    @test isfile("opm.nullmodel.txt")
    @test isfile("opm.scoretest.txt")
    scorepvals = CSV.read("opm.scoretest.txt")[5][1:5]
    @test all(scorepvals .≈ [1.0, 6.45042592e-3, 1.57850425e-5, 4.97907512e-6, 4.56657402e-3])
    rm("opm.nullmodel.txt")
    rm("opm.scoretest.txt")
end

@testset "mask" begin
    polrgwas(@formula(trait ~ 0 + sex), covfile, plkfile, colinds=1:5, outfile="first5snps")
    @test isfile("first5snps.nullmodel.txt")
    @test isfile("first5snps.scoretest.txt")
    @test countlines("first5snps.scoretest.txt") == 6
    scorepvals = CSV.read("first5snps.scoretest.txt")[5]
    @test all(scorepvals .≈ [1.0, 3.00123075e-3, 2.51172149e-5, 1.137311225e-5, 8.31735837e-3])
    rm("first5snps.nullmodel.txt")
    rm("first5snps.scoretest.txt")
end

@testset "test formula" begin
    # score test
    polrgwas(@formula(trait ~ 0 + sex), covfile, plkfile, outfile="GxE", testformula=@formula(trait ~ 0 + snp + snp & sex))
    @test isfile("GxE.nullmodel.txt")
    @test isfile("GxE.scoretest.txt")
    scorepvals = CSV.read("GxE.scoretest.txt")[5][1:5]
    @test all(scorepvals .≈ [1.0, 1.17883709e-2, 1.38956191e-4, 4.44624250e-5, 2.94881322e-2])
    rm("GxE.nullmodel.txt")
    rm("GxE.scoretest.txt")
    # LRT, only first 5 SNPs
    polrgwas(@formula(trait ~ 0 + sex), covfile, plkfile, outfile="GxE", 
    testformula=@formula(trait ~ 0 + snp + snp & sex), test=:LRT, colinds=1:5)
    @test isfile("GxE.nullmodel.txt")
    @test isfile("GxE.lrttest.txt")
    lrtpvals = CSV.read("GxE.lrttest.txt")[end]
    @test all(lrtpvals .≈ [1.0, 7.22410973e-3, 1.01730983e-4, 1.88174211e-5, 2.88295705e-2])
    rm("GxE.nullmodel.txt")
    rm("GxE.lrttest.txt")
end
