"""
    ordinalgwas(nullformula, covfile, plinkfile)
    ordinalgwas(nullformula, df, plinkfile)
    ordinalgwas(fittednullmodel, plinkfile)
    ordinalgwas(fittednullmodel, bedfile, bimfile, bedn)

# Positional arguments 
- `nullformula::Formula`: formula for the null model.
- `covfile::AbstractString`: covariate file (csv) with one header line. One column 
    should be the ordered categorical phenotype coded as integers starting from 1.
- `df::DataFrame`: DataFrame containing response and regressors for null model.
- `plinkfile::AbstractString`: Plink file name without the bed, fam, or bim 
    extensions. If `plinkfile==nothing`, only null model is fitted.  
- `fittednullmodel::StatsModels.DataFrameRegressionModel`: the fitted null model 
    output from `ordinalgwas(nullformula, covfile)` or `ordinalgwas(nullformula, df)`.
- `bedfile::Union{AbstractString,IOStream}`: path to Plink bed file.
- `bimfile::Union{AbstractString,IOStream}`: path to Plink bim file.
- `bedn::Integer`: number of samples in bed file.

# Keyword arguments
- `nullfile::Union{AbstractString, IOStream}`: output file for the fitted null model; default is 
    `ordinalgwas.null.txt`. 
- `pvalfile::Union{AbstractString, IOStream}`: output file for the gwas p-values; default is 
    `ordinalgwas.pval.txt`. 
- `covtype::Vector{DataType}`: type information for `covfile`. This is useful
    when `CSV.read(covarfile)` has parsing errors.  
- `covrowinds::Union{Nothing,AbstractVector{<:Integer}}`: sample indices for covariate file.  
- `testformula::Formula`: formula for test unit. Default is `@formula(trait ~ 0 + snp)`.
- `test::Symbol`: `:score` (default) or `:lrt`.  
- `link::GLM.Link`: `LogitLink()` (default), `ProbitLink()`, `CauchitLink()`,
    or `CloglogLink()`.
- `snpmodel`: `ADDITIVE_MODEL` (default), `DOMINANT_MODEL`, or `RECESSIVE_MODEL`.
- `snpinds::Union{Nothing,AbstractVector{<:Integer}}`: SNP indices for bed file.
- `bedrowinds::Union{Nothing,AbstractVector{<:Integer}}`: sample indices for bed file.
- `solver`: an optimization solver supported by MathProgBase. Default is 
    `NLoptSolver(algorithm=:LD_SLSQP, maxeval=4000)`. Another common choice is 
    `IpoptSolver(print_level=0)`.
- `verbose::Bool`: default is `false`.
"""
function ordinalgwas(
    # positional arguments
    nullformula::Formula,
    covfile::AbstractString,
    plinkfile::Union{Nothing, AbstractString} = nothing;
    # keyword arguments
    covtype::Union{Nothing, Vector{DataType}} = nothing,
    covrowinds::Union{Nothing, AbstractVector{<:Integer}} = nothing,
    kwargs...
    )
    covdf = SnpArrays.makestream(covfile) do io
        CSV.read(io; types=covtype)
    end
    ordinalgwas(nullformula, covrowinds == nothing ? covdf : covdf[covrowinds, :], 
        plinkfile; kwargs...)
end

function ordinalgwas(
    nullformula::Formula,
    nulldf::DataFrame,
    plinkfile::Union{Nothing, AbstractString} = nothing;
    nullfile::Union{AbstractString, IOStream} = "ordinalgwas.null.txt",
    link::GLM.Link = LogitLink(),
    solver = NLoptSolver(algorithm=:LD_SLSQP, maxeval=4000),
    verbose::Bool = false,
    kwargs...
    )
    # fit and output null model
    nm = polr(nullformula, nulldf, link, solver)
    verbose && show(nm)
    SnpArrays.makestream(nullfile, "w") do io
        show(io, nm)
    end
    plinkfile == nothing && (return nm)
    ordinalgwas(nm, plinkfile; solver=solver, verbose=verbose, kwargs...)
end

function ordinalgwas(
    # positional arguments
    fittednullmodel::StatsModels.DataFrameRegressionModel,
    plinkfile::AbstractString;
    # keyword arguments
    testformula::Formula = @eval(@formula($(fittednullmodel.mf.terms.eterms[1]) ~ snp)),
    test::Symbol = :score,
    pvalfile::Union{AbstractString, IOStream} = "ordinalgwas.pval.txt",
    snpmodel::Union{Val{1}, Val{2}, Val{3}} = ADDITIVE_MODEL,
    snpinds::Union{Nothing, AbstractVector{<:Integer}} = nothing,
    bedrowinds::Union{Nothing, AbstractVector{<:Integer}} = nothing,
    solver = NLoptSolver(algorithm=:LD_SLSQP, maxeval=4000),
    verbose::Bool = false
    )
    # locate plink bed, fam, bim file
    if isfile(plinkfile * ".bed")
        bedfile = plinkfile * ".bed"
    else
        fmt = findfirst(isfile, plinkfile * ".bed." .* SnpArrays.ALLOWED_FORMAT)
        fmt == nothing && throw(ArgumentError("bed file not found"))
        bedfile = plinkfile * ".bed." * SnpArrays.ALLOWED_FORMAT[fmt]
    end
    famfile = replace(bedfile, ".bed" => ".fam")
    isfile(famfile) || throw(ArgumentError("fam file not found"))
    bimfile = replace(bedfile, ".bed" => ".bim")
    isfile(bimfile) || throw(ArgumentError("bim file not found"))
    # selected rows should match nobs in null model
    bedn = SnpArrays.makestream(countlines, famfile)
    if bedrowinds == nothing
        nbedrows = bedn
        rowinds = 1:bedn
    else
        nbedrows = eltype(bedrowinds) == Bool ? count(bedrowinds) : length(bedrowinds)
        rowinds = bedrowinds
    end
    nbedrows == nobs(fittednullmodel) || 
        throw(ArgumentError("number of samples in bedrowinds does not match null model"))
    # create SNP mask vector
    if snpinds == nothing
        snpmask = trues(SnpArrays.makestream(countlines, bimfile))
    elseif eltype(snpinds) == Bool
        snpmask = snpinds
    else
        snpmask = falses(SnpArrays.makestream(countlines, bimfile))
        snpmask[snpinds] .= true
    end
    # validate testing method
    test = Symbol(lowercase(string(test)))
    test == :score || test == :lrt || throw(ArgumentError("unrecognized test $test"))
    # gwas
    ordinalgwas(fittednullmodel, bedfile, bimfile, bedn; 
        testformula = testformula, 
        test = test, 
        pvalfile = pvalfile,
        snpmodel = snpmodel, 
        snpmask = snpmask, 
        rowinds = rowinds, 
        solver = solver, 
        verbose = verbose)
    end

function ordinalgwas(
    fittednullmodel::StatsModels.DataFrameRegressionModel,
    bedfile::Union{AbstractString, IOStream}, # full path and bed file name
    bimfile::Union{AbstractString, IOStream}, # full path and bim file name
    bedn::Integer;           # number of samples in bed file
    testformula::Formula = @eval(@formula($(fittednullmodel.mf.terms.eterms[1]) ~ snp)),
    test::Symbol = :score,
    pvalfile::Union{AbstractString, IOStream} = "ordinalgwas.pval.txt", 
    snpmodel::Union{Val{1}, Val{2}, Val{3}} = ADDITIVE_MODEL,
    snpmask::BitVector = trues(SnpArrays.makestream(countlines, bimfile)),
    rowinds::AbstractVector{<:Integer} = 1:bedn, # row indices for SnpArray
    solver = NLoptSolver(algorithm=:LD_SLSQP, maxeval=4000),
    verbose::Bool = false
    )
    # create SnpArray
    genomat = SnpArrays.SnpArray(bedfile, bedn)
    # extra columns in design matrix to be tested
    testdf = fittednullmodel.mf.df # TODO: not type stable here
    testdf[:snp] = zeros(size(fittednullmodel.mf.df, 1))
    mfalt = ModelFrame(testformula, testdf)
    mfalt.terms.intercept = false # drop intercept
    Z = similar(ModelMatrix(mfalt).m)
    # carry out score or LRT test SNP by SNP
    snponly = testformula.rhs == :snp
    cc = SnpArrays.counts(genomat, dims=1) # column counts of genomat
    if test == :score
        ts = OrdinalMultinomialScoreTest(fittednullmodel.model, Z)
        SnpArrays.makestream(pvalfile, "w") do io
            println(io, "chr,pos,snpid,maf,pval")
            SnpArrays.makestream(bimfile) do bimio
                for j in eachindex(snpmask)
                    row = readline(bimio)
                    snpmask[j] || continue
                    maf = (cc[3, j] + 2cc[4, j]) / 2(cc[1, j] + cc[3, j] + cc[4, j])
                    maf > 0.5 && (maf = 1 - maf)
                    if maf == 0 # mono-allelic
                        pval = 1.0
                    else
                        if snponly
                            copyto!(ts.Z, @view(genomat[rowinds, j]), impute=true, model=snpmodel)
                        else # snp + other terms
                            copyto!(testdf[:snp], @view(genomat[rowinds, j]), impute=true, model=snpmodel)
                            mfalt = ModelFrame(testformula, testdf)
                            mfalt.terms.intercept = false # drop intercept
                            ts.Z[:] = ModelMatrix(mfalt).m
                        end
                        pval = polrtest(ts)
                    end
                    snpj = split(row)
                    println(io, "$(snpj[1]),$(snpj[4]),$(snpj[2]),$maf,$pval")
                end
            end
        end
    elseif test == :lrt
        nulldev = deviance(fittednullmodel.model)
        Xaug = [fittednullmodel.model.X Z]
        q = size(Z, 2)
        γ̂ = Vector{Float64}(undef, q) # effect size for columns being tested
        SnpArrays.makestream(pvalfile, "w") do io
            if snponly
                println(io, "chr,pos,snpid,maf,effect,pval")
            else
                print(io, "chr,pos,snpid,maf,")
                for j in 1:q
                    print(io, "effect$j,")
                end
                println(io, "pval")
            end
            SnpArrays.makestream(bimfile) do bimio
                for j in eachindex(snpmask)
                    row = readline(bimio)
                    snpmask[j] || continue
                    maf = (cc[3, j] + 2cc[4, j]) / 2(cc[1, j] + cc[3, j] + cc[4, j])
                    maf > 0.5 && (maf = 1 - maf)
                    if maf == 0 # mono-allelic
                        fill!(γ̂, 0)
                        pval = 1.0
                    else
                        if snponly
                            copyto!(@view(Xaug[:, fittednullmodel.model.p+1]), @view(genomat[rowinds, j]), impute=true, model=snpmodel)
                        else # snp + other terms
                            copyto!(testdf[:snp], @view(genomat[rowinds, j]), impute=true, model=snpmodel)
                            mfalt = ModelFrame(testformula, testdf)
                            mfalt.terms.intercept = false # drop intercept
                            Xaug[:, fittednullmodel.model.p+1:end] = ModelMatrix(mfalt).m
                        end
                        altmodel = polr(Xaug, fittednullmodel.model.Y, fittednullmodel.model.link, solver, wts = fittednullmodel.model.wts)
                        copyto!(γ̂, 1, altmodel.β, fittednullmodel.model.p + 1, q)
                        pval = ccdf(Chisq(q), nulldev - deviance(altmodel))
                    end
                    snpj = split(row)
                    if snponly
                        println(io, "$(snpj[1]),$(snpj[4]),$(snpj[2]),$maf,$(γ̂[1]),$pval")
                    else
                        print(io, "$(snpj[1]),$(snpj[4]),$(snpj[2]),$maf,")
                        for j in 1:q
                            print(io, "$(γ̂[j]),")
                        end
                        println(io, pval)
                    end
                end
            end
        end
    end
    return fittednullmodel
end
