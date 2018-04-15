"""
    polrgwas()

# Keyword arguments

- `plinkfile::AbstractString`: Plink file names without extension.

- `covarfile::AbstractString`: covariate file. One column should be the
ordered categorical phenotype coded as integers starting from 1.

- `outptfile::AbstractString`: output file prefix. Two output files
`prefix.nullmodel.txt` and `prefix.scoretest.txt` will be written.

- `covartype::Vector{DataType}`: type information for `covarfile`. This is useful
when `CSV.read(covarfile)` has parsing errors.

- `test::Symbol`: default `:score`.

- `link::GLM.Link`: `LogitLink()` (default), `ProbitLink()`, `CauchitLink()`,
or `CloglogLink()`
"""
function polrgwas(formula;
    plinkfile::AbstractString = nothing,
    covarfile::AbstractString = nothing,
    covartype::Vector{DataType} = nothing,
    outptfile::AbstractString = "polrgwas",
    test::Symbol = :score,
    link::GLM.Link = LogitLink(),
    verbose::Bool = true
    )

    # fit null model
    covdf = CSV.read(covarfile; types=covartype)
    nm = polr(formula, covdf, link)
    verbose && show(nm)
    nmout = open(outptfile * ".nullmodel.txt", "w")
    show(nmout, nm)
    close(nmout)
    isempty(plinkfile) && (return nothing)

    # carry out score test
    snpdata = SnpData(plinkfile)
    genomat = snpdata.snpmatrix
    ts = PolrScoreTest(nm.model, zeros(snpdata.people, 1))
    outfile = open(outptfile * ".scoretest.txt", "w")
    println(outfile, "chr,pos,snpid,maf,pval")
    for s in 1:snpdata.snps
        copy!(ts.Z, genomat[:, s:s]; impute = true)
        pval = polrtest(ts)
        println(outfile, "$(snpdata.chromosome[s]),$(snpdata.basepairs[s]),$(snpdata.snpid[s]),$(snpdata.maf[s]),$pval")
    end
    close(outfile)
end
