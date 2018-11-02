"""
    polrgwas(formula, covfile, plkfile)

# Position arguments 
- `formula`: a model formula.
- `covfile::AbstractString`: covariate file with header line. One column should be the
ordered categorical phenotype coded as integers starting from 1.
- `plkfile::AbstractString`: Plink file name without the 'bed`, `fam`, or `bim` extension.

# Keyword arguments
- `outfile::AbstractString`: output file prefix; default is "polrgwas". Two output files
`prefix.nullmodel.txt` and `prefix.scoretest.txt` will be written.  
- `covartype::Vector{DataType}`: type information for `covarfile`. This is useful
when `CSV.read(covarfile)` has parsing errors.  
- `test::Symbol`: `:score` or `:lrt`.
- `link::GLM.Link`: `LogitLink()` (default), `ProbitLink()`, `CauchitLink()`,
or `CloglogLink()`
"""
function polrgwas(
    # position arguments
    formula::Formula,
    covfile::AbstractString,
    plkfile::Union{Nothing,AbstractString} = nothing;
    # keyword arguments
    covtype::Union{Nothing,Vector{DataType}} = nothing,
    kwargs...
    )
    covdf = CSV.read(covfile; types=covtype)
    polrgwas(formula, covdf, plkfile; kwargs...)
end

"""
    polrgwas(formula, df, plkfile)

# Position arguments 
- `formula`: a model formula.
- `df::DataFrame`: DataFrame containing response and covariates.
- `plkfile::AbstractString`: Plink file name without the 'bed`, `fam`, or `bim` extension.

# Keyword arguments
- `outfile::AbstractString`: output file prefix; default is "polrgwas". Two output files
`prefix.nullmodel.txt` and `prefix.scoretest.txt` will be written.  
- `covartype::Vector{DataType}`: type information for `covarfile`. This is useful
when `CSV.read(covarfile)` has parsing errors.  
- `test::Symbol`: `:score` or `:LRT`.
- `link::GLM.Link`: `LogitLink()` (default), `ProbitLink()`, `CauchitLink()`,
or `CloglogLink()`
- `colinds::Union{Nothing, AbstractVector{<:Integer}}`: SNP indices.
- `rowinds::Union{Nothing, AbstractVector{<:Integer}}`: sample indices.
"""
function polrgwas(
    formula::Formula,
    df::DataFrame,
    plkfile::Union{Nothing,AbstractString} = nothing;
    # keyword arguments
    outfile::AbstractString = "polrgwas",
    link::GLM.Link = LogitLink(),
    test::Symbol = :score,
    colinds::Union{Nothing, AbstractVector{<:Integer}} = nothing,
    rowinds::Union{Nothing, AbstractVector{<:Integer}} = nothing,
    verbose::Bool = true
    )
    # fit null model
    nm = polr(formula, df, link)
    verbose && show(nm)
    open(outfile * ".nullmodel.txt", "w") do io
        show(io, nm)
    end
    plkfile == nothing && (return nothing)
    # selected rows should match nobs in null model
    rinds = something(rowinds, 1:countlines(plkfile * ".fam"))
    nrows = eltype(rinds) == Bool ? count(rinds) : length(rinds)
    nrows == nobs(nm) || throw(ArgumentError("number of samples SnpArray does not match null model"))
    # create SNP mask vector
    if colinds == nothing 
        cmask = trues(countlines(plkfile * ".bim"))
    elseif eltype(colinds) == Bool
        cmask = colinds
    else
        cmask = falses(countlines(plkfile * ".bim"))
        cmask[colinds] .= true
    end
    # carry out score test
    genomat = SnpArrays.SnpArray(plkfile * ".bed")
    mafreq = SnpArrays.maf(genomat)
    ts = PolrScoreTest(nm.model, zeros(size(genomat, 1), 1))
    open(outfile * "." * string(test) * ".txt", "w") do io
        println(io, "chr,pos,snpid,maf,pval")
        for (j, row) in enumerate(eachline(plkfile * ".bim"))
            cmask[j] || continue
            copyto!(ts.Z, @view(genomat[rinds, j]), impute = true)
            pval = polrtest(ts)
            snpj = split(row)
            println(io, "$(snpj[1]),$(snpj[4]),$(snpj[2]),$(mafreq[j]),$pval")
        end
    end
    nothing
end
