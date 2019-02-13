"""
    ordinalgwas(nullformula, covfile, plinkfile)
    ordinalgwas(nullformula, df, plinkfile)

# Positional arguments 
- `nullformula::Formula`: formula for the null model.
- `covfile::AbstractString`: covariate file with one header line. One column 
    should be the ordered categorical phenotype coded as integers starting from 1.
- `df::DataFrame`: DataFrame containing response and regressors.
- `plinkfile::AbstractString`: Plink file name without the bed, fam, or bim 
    extensions. If `plinkfile==nothing`, only null model is fitted.

# Keyword arguments
- `outfile::AbstractString`: output file prefix; default is `ordinalgwas`. Two CSV output files
    `prefix.nullmodel.txt` and `prefix.scoretest.txt` (or `prefix.lrttest.txt`) will be written.
- `covtype::Vector{DataType}`: type information for `covfile`. This is useful
    when `CSV.read(covarfile)` has parsing errors.  
- `testformula::Formula`: formula for test unit. Default is `@formula(trait ~ 0 + snp)`.
- `test::Symbol`: `:score` (default) or `:LRT`.  
- `link::GLM.Link`: `LogitLink()` (default), `ProbitLink()`, `CauchitLink()`,
    or `CloglogLink()`.
- `snpmodel`: `ADDITIVE_MODEL` (default), `DOMINANT_MODEL`, or `RECESSIVE_MODEL`.
- `colinds::Union{Nothing,AbstractVector{<:Integer}}`: SNP indices.
- `rowinds::Union{Nothing,AbstractVector{<:Integer}}`: sample indices for bed file.
- `solver`: an optimization solver supported by MathProgBase. Default is 
    `NLoptSolver(algorithm=:LD_SLSQP, maxeval=4000)`. Another common choice is 
    `IpoptSolver(print_level=0)`.
- `verbose::Bool`: default is `false`.
"""
function ordinalgwas(
    # positional arguments
    nullformula::Formula,
    covfile::AbstractString,
    plinkfile::Union{Nothing,AbstractString} = nothing;
    # keyword arguments
    covtype::Union{Nothing,Vector{DataType}} = nothing,
    kwargs...
    )
    covdf = CSV.read(covfile; types=covtype)
    ordinalgwas(nullformula, covdf, plinkfile; kwargs...)
end

function ordinalgwas(
    # positional arguments
    nullformula::Formula,
    df::DataFrame,
    plinkfile::Union{Nothing,AbstractString} = nothing;
    # keyword arguments
    testformula::Formula = @eval(@formula($(nullformula.lhs) ~ snp)),
    outfile::AbstractString = "ordinalgwas",
    link::GLM.Link = LogitLink(),
    test::Symbol = :score,
    snpmodel::Union{Val{1}, Val{2}, Val{3}} = ADDITIVE_MODEL,
    colinds::Union{Nothing,AbstractVector{<:Integer}} = nothing,
    rowinds::Union{Nothing,AbstractVector{<:Integer}} = nothing,
    solver = NLoptSolver(algorithm=:LD_SLSQP, maxeval=4000),
    verbose::Bool = false
    )
    # fit null model
    nm = polr(nullformula, df, link, solver)
    verbose && show(nm)
    open(outfile * ".nullmodel.txt", "w") do io
        show(io, nm)
    end
    plinkfile == nothing && (return nothing)
    # selected rows should match nobs in null model
    rinds = something(rowinds, 1:countlines(plinkfile * ".fam"))
    nrows = eltype(rinds) == Bool ? count(rinds) : length(rinds)
    nrows == nobs(nm) || throw(ArgumentError("number of samples in SnpArray does not match null model"))
    # create SNP mask vector
    if colinds == nothing 
        cmask = trues(countlines(plinkfile * ".bim"))
    elseif eltype(colinds) == Bool
        cmask = colinds
    else
        cmask = falses(countlines(plinkfile * ".bim"))
        cmask[colinds] .= true
    end
    # dataframe for alternative model
    dfalt = df
    dfalt[:snp] = zeros(size(df, 1))
    # extra columns in design matrix to be tested
    mfalt = ModelFrame(testformula, dfalt)
    mfalt.terms.intercept = false # drop intercept
    Z = similar(ModelMatrix(mfalt).m)
    # carry out score or LRT test SNP by SNP
    snponly = testformula.rhs == :snp
    genomat = SnpArrays.SnpArray(plinkfile * ".bed")
    mafreq = SnpArrays.maf(genomat) # TODO: need to calibrate according to rowinds
    if test == :score
        ts = OrdinalMultinomialScoreTest(nm.model, Z)
        open(outfile * ".scoretest.txt", "w") do io
            println(io, "chr,pos,snpid,maf,pval")
            for (j, row) in enumerate(eachline(plinkfile * ".bim"))
                cmask[j] || continue
                if mafreq[j] == 0
                    pval = 1.0
                else
                    if snponly
                        copyto!(ts.Z, @view(genomat[rinds, j]), impute = true, model=snpmodel)
                    else # snp + other terms
                        copyto!(dfalt[:snp], @view(genomat[rinds, j]), impute = true, model=snpmodel)
                        mfalt = ModelFrame(testformula, dfalt)
                        mfalt.terms.intercept = false # drop intercept
                        ts.Z[:] = ModelMatrix(mfalt).m
                    end
                    pval = polrtest(ts)
                end
                snpj = split(row)
                println(io, snpj[1], ",", snpj[4], ",", snpj[2], ",", mafreq[j], ",", pval)
            end
        end
    elseif test == :LRT
        nulldev = deviance(nm.model)
        Xaug = [nm.model.X Z]
        q = size(Z, 2)
        γ̂ = Vector{Float64}(undef, q) # effect size for columns being tested
        open(outfile * ".lrttest.txt", "w") do io
            if snponly
                println(io, "chr,pos,snpid,maf,effect,pval")
            else
                print(io, "chr,pos,snpid,maf,")
                for j in 1:q
                    print(io, "effect", j, ",")
                end
                println(io, "pval")
            end
            for (j, row) in enumerate(eachline(plinkfile * ".bim"))
                cmask[j] || continue
                if mafreq[j] == 0
                    fill!(γ̂, 0)
                    pval = 1.0
                else
                    if snponly
                        copyto!(@view(Xaug[:, nm.model.p+1]), @view(genomat[rinds, j]), impute = true, model=snpmodel)
                    else # snp + other terms
                        copyto!(dfalt[:snp], @view(genomat[rinds, j]), impute = true, model=snpmodel)
                        mfalt = ModelFrame(testformula, dfalt)
                        mfalt.terms.intercept = false # drop intercept
                        Xaug[:, nm.model.p+1:end] = ModelMatrix(mfalt).m
                    end
                    altmodel = polr(Xaug, nm.model.Y, nm.model.link, solver, wts = nm.model.wts)
                    copyto!(γ̂, 1, altmodel.β, nm.model.p + 1, q)
                    pval = ccdf(Chisq(q), nulldev - deviance(altmodel))
                end
                snpj = split(row)
                if snponly
                    println(io, snpj[1], ",", snpj[4], ",", snpj[2], ",", mafreq[j], ",", γ̂[1], ",", pval)
                else
                    print(io, snpj[1], ",", snpj[4], ",", snpj[2], ",", mafreq[j], ",")
                    for j in 1:q
                        print(io, γ̂[j], ",")
                    end
                    println(io, pval)
                end
            end
        end
    else
        throw(ArgumentError("unrecognized test $test"))
    end
    nothing
end
