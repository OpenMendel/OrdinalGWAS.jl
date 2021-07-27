"""
    ordinalgwas(nullformula, covfile, geneticfile; kwargs...)
    ordinalgwas(nullformula, df, geneticfile; kwargs...)
    ordinalgwas(fittednullmodel, geneticfile; kwargs...)
    ordinalgwas(fittednullmodel, bedfile, bimfile, bedn; kwargs...)

# Positional arguments 
- `nullformula::FormulaTerm`: formula for the null model.
- `covfile::AbstractString`: covariate file (csv) with one header line. One column 
    should be the ordinal phenotype coded as integers starting from 1.  For example, 
    ordinal phenotypes can be coded as 1, 2, 3, 4 but not 0, 1, 2, 3.  
- `df::DataFrame`: DataFrame containing response and regressors for null model.
- `geneticfile::Union{Nothing, AbstractString}`: File containing genetic information for GWAS.
    This includes a PLINK file name without the .bed, .fam, or .bim 
    extensions or a VCF file without the .vcf extension. If `geneticfile==nothing`, 
    only null model is fitted. If `geneticfile` is provided, bed, bim, and fam file (or vcf) with 
    the same `geneticfile` prefix need to exist. Compressed file formats such as gz and bz2 
    are allowed. Check all allowed formats by `SnpArrays.ALLOWED_FORMAT`. If you're using a VCF file,
    make sure to use the `geneticformat = "VCF"` keyword option, and specificy dosage (:DS) or 
    genotype (:GT) data with the `vcftype` command.
- `fittednullmodel::StatsModels.TableRegressionModel`: the fitted null model 
    output from `ordinalgwas(nullformula, covfile)` or `ordinalgwas(nullformula, df)`.
- `bedfile::Union{AbstractString,IOStream}`: path to Plink bed file with full file name.
- `bimfile::Union{AbstractString,IOStream}`: path to Plink bim file with full file name.
- `bedn::Integer`: number of samples in bed/vcf file.

# Keyword arguments
- `analysistype`::AbstractString: Type of analysis to conduct. Default is `singlesnp`. Other options are `snpset` and `gxe`.
- `geneticformat`::AbstractString: Type of file used for the genetic analysis. `"PLINK"`, `"VCF"`, and `"BGEN"` are currently supported. Default is PLINK.
- `vcftype`::Union{Symbol, Nothing}: Data to extract from the VCF file for the GWAS analysis. `:DS` for dosage or `:GT` for genotypes. Default is nothing.
- `nullfile::Union{AbstractString, IOStream}`: output file for the fitted null model; 
    default is `ordinalgwas.null.txt`. 
- `pvalfile::Union{AbstractString, IOStream}`: output file for the gwas p-values; default is 
    `ordinalgwas.pval.txt`. 
- `covtype::Vector{DataType}`: type information for `covfile`. This is useful
    when `CSV.read(covarfile)` has parsing errors.  
- `covrowinds::Union{Nothing,AbstractVector{<:Integer}}`: sample indices for covariate file.  
- `testformula::FormulaTerm`: formula for test unit. Default is `@formula(trait ~ 0 + snp)`.
- `test::Symbol`: `:score` (default) or `:lrt`.  
- `link::GLM.Link`: `LogitLink()` (default), `ProbitLink()`, `CauchitLink()`,
    or `CloglogLink()`.
- `snpmodel`: `ADDITIVE_MODEL` (default), `DOMINANT_MODEL`, or `RECESSIVE_MODEL`.
- `snpinds::Union{Nothing,AbstractVector{<:Integer}}`: SNP indices for bed/vcf file.
- `geneticrowinds::Union{Nothing,AbstractVector{<:Integer}}`: sample indices for bed/vcf file.
- `solver`: an optimization solver supported by MathProgBase. Default is 
    `NLoptSolver(algorithm=:LD_SLSQP, maxeval=4000)`. Another common choice is 
    `IpoptSolver(print_level=0)`.
- `verbose::Bool`: default is `false`.
- `snpset::Union{Nothing, Integer, AbstractString, AbstractVector{<:Integer}}`: Only include 
    if you are conducting a snpset analysis. An integer indicates a window of SNPs 
    (i.e. every 500 snps). An abstract string allows you to specify an input file, 
    with no header and two columns separated by a space. The first column must contain the snpset ID
    and the second column must contain the snpid's identical to the bimfile. An AbstractVector
    allows you to specify the snps you want to perform one joint snpset test for.
- `e::Union{AbstractString,Symbol}`: Only include if you are conducting a GxE analysis. 
    Enviromental variable to be used to test the GxE interaction.
    For instance, for testing `sex & snp` interaction, use `:sex` or `"sex"`.


# Examples 
The following is an example of basic GWAS with PLINK files:
```julia
plinkfile = "plinkexample"
covfile = "covexample"
ordinalgwas(@formula(trait ~ sex), covfile, plkfile)
```

The following is an example of basic GWAS with a VCF file using dosages then genotypes:
```julia
vcffile = "vcfexample"
covfile = "covexample"
ordinalgwas(@formula(trait ~ sex), covfile, vcffile; 
    geneticfile = "VCF", vcftype = :DS)

ordinalgwas(@formula(trait ~ sex), covfile, vcffile; 
    geneticfile = "VCF", vcftype = :GT)
```

The following is an example of snpset GWAS (every 50 snps). For more types of snpset analyses see documentation:
```julia
ordinalgwas(@formula(trait ~ sex), covfile, plkfile; 
    analysistype = "snpset", snpset = 50)
```

The following is an example of GxE GWAS testing the interaction effect:
```julia
ordinalgwas(@formula(trait ~ sex), covfile, plkfile;
    analysistype = "gxe", e = :sex)
```

"""
function ordinalgwas(
    # positional arguments
    nullformula::FormulaTerm,
    covfile::AbstractString,
    geneticfile::Union{Nothing, AbstractString} = nothing;
    # keyword arguments
    covtype::Union{Nothing, Vector{DataType}} = nothing,
    covrowinds::Union{Nothing, AbstractVector{<:Integer}} = nothing,
    kwargs...
    )
    covdf = SnpArrays.makestream(covfile) do io
        CSV.read(io, DataFrame; types=covtype)
    end
    ordinalgwas(nullformula, covrowinds == nothing ? covdf : covdf[covrowinds, :], 
        geneticfile; kwargs...)
end

function ordinalgwas(
    nullformula::FormulaTerm,
    nulldf::DataFrame,
    geneticfile::Union{Nothing, AbstractString} = nothing;
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
    geneticfile == nothing && (return nm)
    ordinalgwas(nm, geneticfile; solver=solver, verbose=verbose, kwargs...)
end

function ordinalgwas(
    # positional arguments
    fittednullmodel::StatsModels.TableRegressionModel,
    geneticfile::AbstractString;
    # keyword arguments
    analysistype::AbstractString = "singlesnp",
    geneticformat::AbstractString = "PLINK",
    vcftype::Union{Symbol, Nothing} = nothing,
    testformula::FormulaTerm = fittednullmodel.mf.f.lhs ~ Term(:snp),
    test::Symbol = :score,
    pvalfile::Union{AbstractString, IOStream} = "ordinalgwas.pval.txt",
    snpmodel::Union{Val{1}, Val{2}, Val{3}} = ADDITIVE_MODEL,
    snpinds::Union{Nothing, AbstractVector{<:Integer}} = nothing,
    geneticrowinds::Union{Nothing, AbstractVector{<:Integer}} = nothing,
    solver = NLoptSolver(algorithm=:LD_SLSQP, maxeval=4000),
    verbose::Bool = false,
    snpset::Union{Nothing, Integer, AbstractString, #for snpset analysis
        AbstractVector{<:Integer}} = nothing,
    e::Union{Nothing, AbstractString, Symbol} = nothing # for GxE analysis
    )

    # locate plink bed, fam, bim files or VCF file
    lowercase(geneticformat) in ["plink", "vcf", "bgen"] || error("`geneticformat` $geneticformat not valid. Please use 'VCF' or 'PLINK'.")
    isplink = "plink" == lowercase(geneticformat)
    if lowercase(geneticformat) == "plink"
        if isfile(geneticfile * ".bed")
            bedfile = geneticfile * ".bed"
        else
            fmt = findfirst(isfile, geneticfile * ".bed." .* SnpArrays.ALLOWED_FORMAT)
            fmt == nothing && throw(ArgumentError("bed file not found"))
            bedfile = geneticfile * ".bed." * SnpArrays.ALLOWED_FORMAT[fmt]
        end
        famfile = replace(bedfile, ".bed" => ".fam")
        isfile(famfile) || throw(ArgumentError("fam file not found"))
        bimfile = replace(bedfile, ".bed" => ".bim")
        isfile(bimfile) || throw(ArgumentError("bim file not found"))
        # selected rows should match nobs in null model
        bedn = SnpArrays.makestream(countlines, famfile)
    elseif lowercase(geneticformat) == "vcf"
        vcftype in [:GT, :DS] || throw(ArgumentError("vcftype not specified. Allowable types are :GT for genotypes and :DS for dosages."))
        if isfile(geneticfile * ".vcf")
            vcffile = geneticfile * ".vcf"
        else
            fmt = findfirst(isfile, geneticfile * ".vcf." .* SnpArrays.ALLOWED_FORMAT)
            fmt == nothing && throw(ArgumentError("VCF file not found"))
            vcffile = geneticfile * ".vcf." * SnpArrays.ALLOWED_FORMAT[fmt]
        end
        bedn = VCFTools.nsamples(vcffile)
    elseif lowercase(geneticformat) == "bgen"
        if isfile(geneticfile * ".bgen")
            bgenfile = geneticfile * ".bgen"
        else
            fmt = findfirst(isfile, geneticfile * ".bgen." .* SnpArrays.ALLOWED_FORMAT)
            fmt == nothing && throw(ArgumentError("BGEN file not found"))
            bgenfile = geneticfile * ".bgen." * SnpArrays.ALLOWED_FORMAT[fmt]
        end
        b = Bgen(bgenfile)
        bedn = n_samples(b)
    end
    if geneticrowinds == nothing
        nbedrows = bedn
        rowinds = 1:bedn
    else
        nbedrows = eltype(geneticrowinds) == Bool ? count(geneticrowinds) : length(geneticrowinds)
        rowinds = geneticrowinds
    end

    nbedrows == nobs(fittednullmodel) || 
        throw(ArgumentError("number of samples in geneticrowinds does not match null model"))

    # validate testing method
    test = Symbol(lowercase(string(test)))
    test == :score || test == :lrt || throw(ArgumentError("unrecognized test $test"))

    # gwas
    if lowercase(geneticformat) == "plink" #plink
        ordinalgwas(fittednullmodel, bedfile, bimfile, bedn;
            analysistype = analysistype,
            testformula = testformula, 
            test = test, 
            pvalfile = pvalfile,
            snpmodel = snpmodel, 
            snpinds = snpinds, 
            bedrowinds = rowinds, 
            solver = solver, 
            verbose = verbose,
            snpset = snpset,
            e = e)
    elseif lowercase(geneticformat) == "vcf" #vcf
        ordinalgwas(fittednullmodel, vcffile, bedn, vcftype; 
            analysistype = analysistype,
            testformula = testformula, 
            test = test, 
            pvalfile = pvalfile,
            snpmodel = snpmodel, 
            snpinds = snpinds, 
            vcfrowinds = rowinds, 
            solver = solver, 
            verbose = verbose,
            snpset = snpset,
            e = e)
    else #bgen
        ordinalgwas(fittednullmodel, bgenfile, bedn; 
            analysistype = analysistype,
            testformula = testformula, 
            test = test, 
            pvalfile = pvalfile,
            snpmodel = snpmodel, 
            snpinds = snpinds, 
            bgenrowinds = rowinds, 
            solver = solver, 
            verbose = verbose,
            snpset = snpset,
            e = e)
    end
end

# For PLINK Analysis
function ordinalgwas(
    fittednullmodel::StatsModels.TableRegressionModel,
    bedfile::Union{AbstractString, IOStream}, # full path and bed file name
    bimfile::Union{AbstractString, IOStream}, # full path and bim file name
    bedn::Integer;           # number of samples in bed file
    analysistype::AbstractString = "singlesnp",
    testformula::FormulaTerm = fittednullmodel.mf.f.lhs ~ Term(:snp),
    test::Symbol = :score,
    pvalfile::Union{AbstractString, IOStream} = "ordinalgwas.pval.txt", 
    snpmodel::Union{Val{1}, Val{2}, Val{3}} = ADDITIVE_MODEL,
    snpinds::Union{Nothing, AbstractVector{<:Integer}} = nothing,
    bedrowinds::AbstractVector{<:Integer} = 1:bedn, # row indices for SnpArray
    solver = NLoptSolver(algorithm=:LD_SLSQP, maxeval=4000),
    verbose::Bool = false,
    snpset::Union{Nothing, Integer, AbstractString, #for snpset analysis
        AbstractVector{<:Integer}} = nothing,
    e::Union{Nothing, AbstractString, Symbol} = nothing # for GxE analysis
    )
    # create SnpArray
    genomat = SnpArrays.SnpArray(bedfile, bedn)
    cc = SnpArrays.counts(genomat, dims=1) # column counts of genomat
    mafs = SnpArrays.maf(genomat)

    # create SNP mask vector
    if snpinds == nothing
        snpmask = trues(SnpArrays.makestream(countlines, bimfile))
    elseif eltype(snpinds) == Bool
        snpmask = snpinds
    else
        snpmask = falses(SnpArrays.makestream(countlines, bimfile))
        snpmask[snpinds] .= true
    end

    analysistype = lowercase(analysistype)
    analysistype in ["singlesnp", "snpset", "gxe"] || error("Analysis type $analysis invalid option. 
    Available options are 'singlesnp', 'snpset' and 'gxe'.")

    # determine analysis type
    if analysistype == "singlesnp"

        # carry out score or LRT test SNP by SNP
        snponly = testformula.rhs == Term(:snp)

        # extra columns in design matrix to be tested
        testdf = DataFrame(fittednullmodel.mf.data) # TODO: not type stable here
        testdf[!, :snp] = zeros(size(fittednullmodel.mm, 1))
        Z = similar(modelmatrix(testformula, testdf))
        SnpArrays.makestream(pvalfile, "w") do io
            if test == :score
                println(io, "chr,pos,snpid,maf,hwepval,pval")
                ts = OrdinalMultinomialScoreTest(fittednullmodel.model, Z)
            else # lrt 
                nulldev = deviance(fittednullmodel.model)
                Xaug = [fittednullmodel.model.X Z]
                q = size(Z, 2)
                γ̂ = Vector{Float64}(undef, q) # effect size for columns being tested
                if snponly
                    println(io, "chr,pos,snpid,maf,hwepval,effect,stder,pval")
                else
                    print(io, "chr,pos,snpid,maf,hwepval,")
                    for j in 1:q
                        print(io, "effect$j,")
                    end
                    println(io, "pval")
                end
            end
            SnpArrays.makestream(bimfile) do bimio
                for j in eachindex(snpmask)
                    row = readline(bimio)
                    snpmask[j] || continue
                    snpj = split(row)
                    hwepval = SnpArrays.hwe(cc[1, j], cc[3, j], cc[4, j])
                    maf = mafs[j]
                    if test == :score
                        if maf == 0 # mono-allelic
                            pval = 1.0
                        else
                            if snponly
                                copyto!(ts.Z, @view(genomat[bedrowinds, j]), 
                                impute=true, model=snpmodel)
                            else # snp + other terms
                                copyto!(testdf[!, :snp], @view(genomat[bedrowinds, j]),
                                    impute = true, model = snpmodel)
                                ts.Z[:] = modelmatrix(testformula, testdf)
                            end
                            pval = polrtest(ts)
                        end
                        println(io, "$(snpj[1]),$(snpj[4]),$(snpj[2]),",
                            "$maf,$hwepval,$pval")
                    elseif test == :lrt
                        if maf == 0 # mono-allelic
                            fill!(γ̂, 0)
                            stderr = -1.0
                            pval = 1.0
                        else
                            if snponly
                                copyto!(@view(Xaug[:, fittednullmodel.model.p+1]), 
                                    @view(genomat[bedrowinds, j]), 
                                    impute=true, model=snpmodel)
                            else # snp + other terms
                                copyto!(testdf[!, :snp], @view(genomat[bedrowinds, j]), 
                                    impute=true, model=snpmodel)
                                Xaug[:, fittednullmodel.model.p+1:end] = modelmatrix(testformula, testdf)
                            end
                            altmodel = polr(Xaug, fittednullmodel.model.Y, 
                                fittednullmodel.model.link, solver, 
                                wts = fittednullmodel.model.wts)
                            copyto!(γ̂, 1, altmodel.β, fittednullmodel.model.p + 1, q)
                            pval = ccdf(Chisq(q), nulldev - deviance(altmodel))
                            stderr = stderror(altmodel)[fittednullmodel.model.npar + 1]
                        end
                        if snponly
                            println(io, "$(snpj[1]),$(snpj[4]),$(snpj[2]),$maf,$hwepval,",
                                "$(γ̂[1]),$stderr,$pval")
                        else
                            print(io, "$(snpj[1]),$(snpj[4]),$(snpj[2]),$maf,$hwepval,")
                            for j in 1:q
                                print(io, "$(γ̂[j]),")
                            end
                            println(io, pval)
                        end
                    end
                end
            end
        end
    elseif analysistype == "snpset"
        # determine type of snpset analysis 
        if isa(snpset, Nothing)
            @warn("Nothing set for `snpset`. 
            This will default to `singlesnp` analysis (windowsize = 1).")
            setlength = 1
        elseif isa(snpset, AbstractString)
            isfile(snpset) || throw(ArgumentError("snpset file not found, 
            to specify a window replace snpset string with a window size"))
            #first column SNPset ID, second column SNP ID
            snpsetFile = CSV.read(snpset, DataFrame, header = [:snpset_id, :snp_id], delim = " ")
            #make sure it matches bim file
            biminfo = CSV.read(bimfile, DataFrame, header = [:chr, :snp_id, :c3, :bp, :c5, :c6], delim = "\t")
            snpsetFile[!, :snp_id] == biminfo[!, :snp_id] || throw(ArgumentError("snp order in snpset file
            must match (in the same order) bimfile")) 
            snpset_ids = unique(snpsetFile[!, :snpset_id])
            nSets = length(snpset_ids)
            setlength = 0
        elseif isa(snpset, Integer)
            setlength = snpset
        else #abstract vector (boolean of true at indicies or range or indicies)
            setlength = -1
        end

        # conduct analysis based on type 
        if setlength > 0 #single snp analysis or window
            Z = zeros(size(fittednullmodel.mm, 1), setlength) # column counts of genomat
            totalsnps = SnpArrays.makestream(countlines, bimfile)
            SnpArrays.makestream(pvalfile, "w") do io
                if test == :score
                    println(io, "startchr,startpos,startsnpid,endchr,endpos,endsnpid,pval")
                    ts = OrdinalMultinomialScoreTest(fittednullmodel.model, Z)
                else # lrt 
                    println(io, "startchr,startpos,startsnpid,endchr,",
                    "endpos,endsnpid,l2normeffect,pval")
                    nulldev = deviance(fittednullmodel.model)
                    Xaug = [fittednullmodel.model.X Z]
                    γ̂ = Vector{Float64}(undef, setlength) # effect size for columns being tested
                end
                SnpArrays.makestream(bimfile) do bimio
                    q = setlength
                    for j in 1:q:totalsnps
                        endj = j + q - 1  
                        rowj = readline(bimio)  
                        if endj >= totalsnps
                            endj = totalsnps
                            q = totalsnps - j + 1
                            #length of Z will be different
                            if test == :score 
                                Z = zeros(size(fittednullmodel.mm, 1), q)
                                ts = OrdinalMultinomialScoreTest(fittednullmodel.model, Z)
                            else
                                Xaug = [fittednullmodel.model.X zeros(size(
                                    fittednullmodel.mm, 1), q)]
                            end
                        end
                        for i in 1:(q - 2) #
                            readline(bimio)
                        end
                        endj == totalsnps ? rowj_s = rowj : rowj_s = readline(bimio)
                        snpj = split(rowj)
                        snpj_s = split(rowj_s)
                        if test == :score 
                            if all(@view(mafs[j:endj]) .== 0) # all mono-allelic, unlikely but just in case
                                pval = 1.0
                            else
                                copyto!(ts.Z, @view(genomat[bedrowinds, j:endj]), impute=true, model=snpmodel)
                                pval = polrtest(ts)
                            end
                            println(io, "$(snpj[1]),$(snpj[4]),$(snpj[2]),$(snpj_s[1]),",
                                "$(snpj_s[4]),$(snpj_s[2]),$pval")
                        elseif test == :lrt 
                            if all(@view(mafs[j:endj]) .== 0) # all mono-allelic, unlikely but just in case
                                fill!(γ̂, 0)
                                pval = 1.0
                            else
                                copyto!(@view(Xaug[:, (fittednullmodel.model.p+1):end]), 
                                        @view(genomat[bedrowinds, j:endj]), 
                                        impute=true, model=snpmodel)
                                altmodel = polr(Xaug, fittednullmodel.model.Y, 
                                    fittednullmodel.model.link, solver, 
                                    wts = fittednullmodel.model.wts)
                                copyto!(γ̂, @view(altmodel.β[(fittednullmodel.model.p+1):end]))#, fittednullmodel.model.p + 1, setlength)
                                l2normeffect = norm(γ̂)
                                pval = ccdf(Chisq(q), nulldev - deviance(altmodel))
                            end
                            println(io, "$(snpj[1]),$(snpj[4]),$(snpj[2]),",
                                "$(snpj_s[1]),$(snpj_s[4]),$(snpj_s[2]),$l2normeffect,$pval")
                        end
                    end
                end
            end
        elseif setlength == 0 #snpset is defined by snpset file
            SnpArrays.makestream(pvalfile, "w") do io
                test == :score ? println(io, "snpsetid,nsnps,pval") : println(io, 
                    "snpsetid,nsnps,l2normeffect,pval")
                for j in eachindex(snpset_ids)
                    snpset_id = snpset_ids[j]
                    snpinds = findall(snpsetFile[!, :snpset_id] .== snpset_id)
                    q = length(snpinds)
                    Z = zeros(size(fittednullmodel.mm, 1), q)
                    γ̂ = Vector{Float64}(undef, q)
                    Xaug = [fittednullmodel.model.X Z]
                    if all(@view(mafs[snpinds]) .== 0) # all mono-allelic, unlikely but just in case
                        pval = 1.0
                        l2normeffect = 0.0
                        test == :score ? println(io, "$(snpset_id),$q,$pval") : 
                        println(io, "$(snpset_id),$q,$l2normeffect,$pval")
                    elseif test == :score
                        ts = OrdinalMultinomialScoreTest(fittednullmodel.model, Z)
                        copyto!(ts.Z, @view(genomat[bedrowinds, snpinds]), impute=true,
                            model=snpmodel)
                        pval = polrtest(ts)
                        println(io, "$(snpset_id),$q,$pval")
                    elseif test == :lrt
                        nulldev = deviance(fittednullmodel.model)
                        copyto!(@view(Xaug[:, fittednullmodel.model.p+1:end]), 
                                @view(genomat[bedrowinds, snpinds]), 
                                impute=true, model=snpmodel)
                        altmodel = polr(Xaug, fittednullmodel.model.Y, 
                            fittednullmodel.model.link, solver, 
                            wts = fittednullmodel.model.wts)
                        copyto!(γ̂, 1, altmodel.β, fittednullmodel.model.p + 1, q)
                        l2normeffect = norm(γ̂)
                        pval = ccdf(Chisq(q), nulldev - deviance(altmodel))
                        println(io, "$(snpset_id),$q,$l2normeffect,$pval")
                    end
                end
            end
        else #setlength == -1 (testing just one set with specified snps in snpset)
            snpset = eltype(snpset) == Bool ? findall(snpset) : snpset 
            SnpArrays.makestream(pvalfile, "w") do io
                if all(@view(mafs[snpset]) .== 0) # all mono-allelic, unlikely but just in case
                    pval = 1.0
                    l2normeffect = 0.0
                    test == :score ? println(io, "The joint pvalue of snps indexed",
                    " at $(snpset) is $pval") : println(io, "The l2norm of the effect size vector",
                    " is $l2normeffect and joint pvalue of snps indexed", 
                    " at $(snpset) is $pval")
                else
                    q = length(snpset)
                    γ̂ = Vector{Float64}(undef, q)
                    Z = zeros(size(fittednullmodel.mm, 1), q)
                    if test == :score
                        ts = OrdinalMultinomialScoreTest(fittednullmodel.model, Z)
                        copyto!(ts.Z, @view(genomat[bedrowinds, snpset]), impute=true, model=snpmodel)
                        pval = polrtest(ts)
                        println(io, "The joint pvalue of snps indexed",
                         " at $(snpset) is $pval")
                    elseif test == :lrt
                        nulldev = deviance(fittednullmodel.model)
                        Xaug = [fittednullmodel.model.X Z]
                        copyto!(@view(Xaug[:, fittednullmodel.model.p+1:end]), 
                                @view(genomat[bedrowinds, snpset]), 
                                impute=true, model=snpmodel)
                        altmodel = polr(Xaug, fittednullmodel.model.Y, 
                            fittednullmodel.model.link, solver, 
                            wts = fittednullmodel.model.wts)
                        copyto!(γ̂, 1, altmodel.β, fittednullmodel.model.p + 1, q)
                        l2normeffect = norm(γ̂)
                        pval = ccdf(Chisq(q), nulldev - deviance(altmodel))
                        println(io, "The l2norm of the effect size vector",
                        " is $l2normeffect and joint pvalue of snps indexed", 
                        " at $(snpset) is $pval")
                    end
                end
            end
        end
    else #analysistype == "gxe"
        isnothing(e) && 
            @error("GxE analysis indicated but not environmental variable keyword argument: `e` set.")
        Xaug = [fittednullmodel.model.X zeros(size(fittednullmodel.mm, 1))]
        Xaug2 = [fittednullmodel.model.X zeros(size(fittednullmodel.mm, 1), 2)]
        envvar = modelmatrix(FormulaTerm(fittednullmodel.mf.f.lhs, Term(Symbol(e))),
                 DataFrame(fittednullmodel.mf.data))
        testvec = Matrix{Float64}(undef, size(envvar))
        snpeffectnull = 0.0
        SnpArrays.makestream(pvalfile, "w") do io
            if test == :score 
                println(io, "chr,pos,snpid,maf,hwepval,snpeffectnull,pval")
            else 
                γ̂ = 0.0 # effect size for columns being tested
                println(io, "chr,pos,snpid,maf,hwepval,snpeffectnull,snpeffectfull,GxEeffect,pval")
            end
            SnpArrays.makestream(bimfile) do bimio
                for j in eachindex(snpmask)
                    row = readline(bimio)
                    snpmask[j] || continue
                    hwepval = SnpArrays.hwe(cc[1, j], cc[3, j], cc[4, j])
                    maf = mafs[j]
                    snpj = split(row)
                    if maf == 0 # mono-allelic
                        γ̂ = 0.0
                        pval = 1.0
                        snpeffectfull = 0.0
                        snpeffectnull = 0.0
                    elseif test == :score
                        copyto!(@view(Xaug[:, end]), @view(genomat[bedrowinds,
                        j]), impute=true, model=snpmodel)
                        copyto!(testvec, @view(Xaug[:, end]) .* envvar)
                        nm = polr(Xaug, fittednullmodel.model.Y, 
                            fittednullmodel.model.link, solver, wts = fittednullmodel.model.wts)
                        if Inf in nm.vcov #singular design matrix
                            snpeffectnull = 0.0
                            pval = 1.0
                        else
                            snpeffectnull = nm.β[end]
                            ts = OrdinalMultinomialScoreTest(nm, testvec)
                            pval = polrtest(ts)
                        end
                    elseif test == :lrt 
                        copyto!(@view(Xaug[:, end]), @view(genomat[bedrowinds,
                            j]), impute=true, model=snpmodel)
                        copyto!(@view(Xaug2[:, end - 1]), @view(Xaug[:, end]))
                        copyto!(@view(Xaug2[:, end]), @view(Xaug[:, end]) .*
                            envvar)
                        nm = polr(Xaug, fittednullmodel.model.Y, 
                        fittednullmodel.model.link, solver, wts = fittednullmodel.model.wts)
                        if Inf in nm.vcov #singular design matrix
                            snpeffectnull = 0.0
                            snpeffectfull = 0.0
                            γ̂ = 0.0
                            pval = 1.0
                        else
                            snpeffectnull = nm.β[end]
                            nulldev = deviance(nm)
                            altmodel = polr(Xaug2, fittednullmodel.model.Y, 
                                fittednullmodel.model.link, solver, 
                                wts = fittednullmodel.model.wts)
                            γ̂ = altmodel.β[end]
                            snpeffectfull = altmodel.β[end-1]
                            pval = ccdf(Chisq(1), nulldev - deviance(altmodel))
                        end
                    end
                    if test == :score
                        println(io, "$(snpj[1]),$(snpj[4]),$(snpj[2]),$maf,",
                        "$hwepval,$snpeffectnull,$pval")
                    else
                        println(io, "$(snpj[1]),$(snpj[4]),$(snpj[2]),$maf,",
                        "$hwepval,$snpeffectnull,$snpeffectfull,$γ̂,$pval")
                    end
                end
            end
        end
    end
    return fittednullmodel
end


# For VCF Analysis
function ordinalgwas(
    fittednullmodel::StatsModels.TableRegressionModel,
    vcffile::Union{AbstractString, IOStream}, # full path and vcf file name
    nsamples::Integer,          # number of samples in bed file
    vcftype::Symbol;           # :GT = genotype, :DS = dosage
    analysistype::AbstractString = "singlesnp",
    testformula::FormulaTerm = fittednullmodel.mf.f.lhs ~ Term(:snp),
    test::Symbol = :score,
    pvalfile::Union{AbstractString, IOStream} = "ordinalgwas.pval.txt", 
    snpmodel::Union{Val{1}, Val{2}, Val{3}} = ADDITIVE_MODEL,
    snpinds::Union{Nothing, AbstractVector{<:Integer}} = nothing,
    vcfrowinds::AbstractVector{<:Integer} = 1:nsamples, # row indices for VCF array
    solver = NLoptSolver(algorithm=:LD_SLSQP, maxeval=4000),
    verbose::Bool = false,
    snpset::Union{Nothing, Integer, AbstractString, #for snpset analysis
        AbstractVector{<:Integer}} = nothing,
    e::Union{Nothing, AbstractString, Symbol} = nothing # for GxE analysis
    )

    # get number of SNPs in file
    nsnps = nrecords(vcffile)

    # for VCFTools, snpmodel is coded differently 
    snpmodel = modelingdict[snpmodel]

    # create SNP mask vector
    if snpinds == nothing
        snpmask = trues(nsnps)
    elseif eltype(snpinds) == Bool
        snpmask = snpinds
    else
        snpmask = falses(nsnps)
        snpmask[snpinds] .= true
    end

    analysistype = lowercase(analysistype)
    analysistype in ["singlesnp", "snpset", "gxe"] || error("Analysis type $analysis invalid option. 
    Available options are 'singlesnp', 'snpset' and 'gxe'.")
    # open VCF File 
    reader = VCF.Reader(openvcf(vcffile))

    # determine analysis type
    if analysistype == "singlesnp"
        # extra columns in design matrix to be tested
        testdf = DataFrame(fittednullmodel.mf.data) # TODO: not type stable here
        testdf[!, :snp] = zeros(size(fittednullmodel.mm, 1))
        Z = similar(modelmatrix(testformula, testdf))

        # create holders for chromome, position, id, dosage/snps 
        rec_chr = Array{Any, 1}(undef, 1)
        rec_pos = Array{Any, 1}(undef, 1)
        rec_ids = Array{Any, 1}(undef, 1)
        gholder = zeros(Union{Missing, Float64}, nsamples)

        # carry out score or LRT test SNP by SNP
        snponly = testformula.rhs == Term(:snp)
        SnpArrays.makestream(pvalfile, "w") do io
            if test == :score 
                println(io, "chr,pos,snpid,pval")
                ts = OrdinalMultinomialScoreTest(fittednullmodel.model, Z)
            else 
                nulldev = deviance(fittednullmodel.model)
                Xaug = [fittednullmodel.model.X Z]
                q = size(Z, 2)
                γ̂ = Vector{Float64}(undef, q) # effect size for columns being tested
                if snponly
                    println(io, "chr,pos,snpid,effect,stder,pval")
                else
                    print(io, "chr,pos,snpid,")
                    for j in 1:q
                        print(io, "effect$j,")
                    end
                    println(io, "pval")
                end        
            end
            for j in eachindex(snpmask)
                if !snpmask[j] #skip snp, must read marker still. 
                    copy_gt!(gholder, reader; model = snpmodel, impute = true,
                    record_chr = rec_chr, record_pos = rec_pos, record_ids = rec_ids)
                    continue
                end
                if vcftype == :GT #genotype 
                    copy_gt!(gholder, reader; model = snpmodel, impute = true,
                    record_chr = rec_chr, record_pos = rec_pos, record_ids = rec_ids)
                else #dosage
                    copy_ds!(gholder, reader; model = snpmodel, impute = true,
                    record_chr = rec_chr, record_pos = rec_pos, record_ids = rec_ids)
                end
                if test == :score
                    if snponly
                        copyto!(ts.Z, @view(gholder[vcfrowinds]))
                    else # snp + other terms
                        copyto!(testdf[!, :snp], @view(gholder[vcfrowinds]))
                        ts.Z[:] = modelmatrix(testformula, testdf)
                    end
                    pval = polrtest(ts)
                    println(io, "$(rec_chr[1]),$(rec_pos[1]),$(rec_ids[1][1]),$pval")
                elseif test == :lrt 
                    if snponly
                        copyto!(@view(Xaug[:, fittednullmodel.model.p+1]), 
                            @view(gholder[vcfrowinds]))
                    else # snp + other terms
                        copyto!(testdf[!, :snp], @view(gholder[vcfrowinds]))
                        Xaug[:, fittednullmodel.model.p+1:end] = modelmatrix(testformula, testdf)
                    end
                    altmodel = polr(Xaug, fittednullmodel.model.Y, 
                        fittednullmodel.model.link, solver, 
                        wts = fittednullmodel.model.wts)
                    copyto!(γ̂, 1, altmodel.β, fittednullmodel.model.p + 1, q)
                    pval = ccdf(Chisq(q), nulldev - deviance(altmodel))
                    stderr = stderror(altmodel)[fittednullmodel.model.npar + 1]
                    if snponly
                        println(io, "$(rec_chr[1]),$(rec_pos[1]),$(rec_ids[1][1]),$(γ̂[1]),$stderr,$pval")
                    else
                        print(io, "$(rec_chr[1]),$(rec_pos[1]),$(rec_ids[1][1]),")
                        for j in 1:q
                            print(io, "$(γ̂[j]),")
                        end
                        println(io, pval)
                    end
                end
            end
        end
    elseif analysistype == "snpset"
        # max size of a snpset length
        maxsnpset = 1

        #determine snpset
        if isa(snpset, Nothing)
            setlength = 1
            maxsnpset = 1
        elseif isa(snpset, AbstractString)
            isfile(snpset) || throw(ArgumentError("snpset file not found, 
            to specify a window replace snpset string with a window size"))
            #first column SNPset ID, second column SNP ID
            snpsetFile = CSV.read(snpset, DataFrame, header = [:snpset_id, :snp_id], delim = " ")
            maxsnpset = combine(groupby(snpsetFile, :snpset_id), :snp_id => length => :snpset_length) |> 
                x -> maximum(x.snpset_length)
            snpset_ids = unique(snpsetFile[!, :snpset_id])
            nSets = length(snpset_ids)
            setlength = 0
        elseif isa(snpset, Integer)
            setlength = snpset
            maxsnpset = snpset 
        else #abstract vector (boolean of true at indicies or range or indicies)
            setlength = -1
            maxsnpset = count(snpset .!= 0)
        end

        # create holders for chromome, position, id, dosage/gt
        rec_chr = Array{Any, 1}(undef, maxsnpset)
        rec_pos = Array{Any, 1}(undef, maxsnpset)
        rec_ids = Array{Any, 1}(undef, maxsnpset)
        gholder = zeros(Union{Missing, Float64}, nsamples, maxsnpset)

        if setlength > 0 #single snp analysis or window
            Z = zeros(size(fittednullmodel.mm, 1), setlength) #
            q = setlength
            SnpArrays.makestream(pvalfile, "w") do io
                if test == :score
                    println(io, "startchr,startpos,startsnpid,endchr,endpos,endsnpid,pval")
                    ts = OrdinalMultinomialScoreTest(fittednullmodel.model, Z)
                elseif test == :lrt 
                    println(io, "startchr,startpos,startsnpid,endchr,",
                    "endpos,endsnpid,l2normeffect,pval")
                    nulldev = deviance(fittednullmodel.model)
                    Xaug = [fittednullmodel.model.X Z]
                    γ̂ = Vector{Float64}(undef, setlength) # effect size for columns being tested
                end
                for j in 1:q:nsnps
                    endj = j + q - 1    
                    if endj >= nsnps
                        endj = nsnps
                        q = nsnps - j + 1
                        #length of Z will be different
                        gholder = zeros(Union{Missing, Float64}, nsamples, q)
                        rec_chr = Array{Any, 1}(undef, q)
                        rec_pos = Array{Any, 1}(undef, q)
                        rec_ids = Array{Any, 1}(undef, q)
                        if test == :score
                            Z = zeros(size(fittednullmodel.mm, 1), q)
                            ts = OrdinalMultinomialScoreTest(fittednullmodel.model, Z)
                        elseif test == :lrt 
                            Xaug = [fittednullmodel.model.X zeros(size(
                            fittednullmodel.mm, 1), q)]
                        end
                    end
                    if vcftype == :GT #genotype 
                        copy_gt!(gholder, reader; model = snpmodel, impute = true,
                        record_chr = rec_chr, record_pos = rec_pos, record_ids = rec_ids)
                    else #dosage
                        copy_ds!(gholder, reader; model = snpmodel, impute = true,
                        record_chr = rec_chr, record_pos = rec_pos, record_ids = rec_ids)
                    end
                    if test == :score
                        copyto!(ts.Z, @view(gholder[vcfrowinds, :]))
                        pval = polrtest(ts)
                        println(io, "$(rec_chr[1]),$(rec_pos[1]),$(rec_ids[1][1]),",
                        "$(rec_chr[end]),$(rec_pos[end]),$(rec_ids[end][end]),$pval")
                    elseif test == :lrt 
                        copyto!(@view(Xaug[:, (fittednullmodel.model.p+1):end]), 
                        @view(gholder[vcfrowinds]))
                        altmodel = polr(Xaug, fittednullmodel.model.Y, 
                        fittednullmodel.model.link, solver, 
                        wts = fittednullmodel.model.wts)
                        copyto!(γ̂, @view(altmodel.β[(fittednullmodel.model.p+1):end]))#, fittednullmodel.model.p + 1, setlength)
                        l2normeffect = norm(γ̂)
                        pval = ccdf(Chisq(q), nulldev - deviance(altmodel))
                        println(io, "$(rec_chr[1]),$(rec_pos[1]),$(rec_ids[1][1]),",
                        "$(rec_chr[end]),$(rec_pos[end]),$(rec_ids[end][end]),$l2normeffect,$pval")
                    end
                end
            end
        elseif setlength == 0 #snpset is defined by snpset file
            @warn("This method requires reading in the entire VCF File.
             This can take a lot of memory for large files, as they must be brought into memory.")
            if vcftype == :GT #genotype 
                genomat = convert_gt(Float64, vcffile; 
                model = snpmodel, impute = true, 
                center = false, scale = false)
            else #dosage
                genomat = convert_ds(Float64, vcffile; model = snpmodel,
                key="DS", impute = true, center = false, scale = false)
            end
            SnpArrays.makestream(pvalfile, "w") do io
                test == :score ? println(io, "snpsetid,nsnps,pval") : println(io, 
                    "snpsetid,nsnps,l2normeffect,pval")
                for j in eachindex(snpset_ids)
                    snpset_id = snpset_ids[j]
                    snpinds = findall(snpsetFile[!, :snpset_id] .== snpset_id)
                    q = length(snpinds)
                    Z = zeros(size(fittednullmodel.mm, 1), q)
                    if test == :score
                        ts = OrdinalMultinomialScoreTest(fittednullmodel.model, Z)
                        copyto!(ts.Z, @view(genomat[vcfrowinds, snpinds]))
                        pval = polrtest(ts)
                        println(io, "$(snpset_id),$q,$pval")
                    elseif test == :lrt
                        γ̂ = Vector{Float64}(undef, q)
                        Xaug = [fittednullmodel.model.X Z]
                        nulldev = deviance(fittednullmodel.model)
                        copyto!(@view(Xaug[:, fittednullmodel.model.p+1:end]), 
                                @view(genomat[vcfrowinds, snpinds]))
                        altmodel = polr(Xaug, fittednullmodel.model.Y, 
                            fittednullmodel.model.link, solver, 
                            wts = fittednullmodel.model.wts)
                        copyto!(γ̂, 1, altmodel.β, fittednullmodel.model.p + 1, q)
                        l2normeffect = norm(γ̂)
                        pval = ccdf(Chisq(q), nulldev - deviance(altmodel))
                        println(io, "$(snpset_id),$q,$l2normeffect,$pval")
                    end
                end
            end
        else #setlength == -1 (testing just one set with specified snps in snpset)
            @warn("This method requires reading in the entire VCF File.
            This can take a lot of memory for large files, as they must be brought into memory.")
            snpset = eltype(snpset) == Bool ? findall(snpset) : snpset 
            if vcftype == :GT #genotype 
                genomat = convert_gt(Float64, vcffile; 
                model = snpmodel, impute = true, 
                center = false, scale = false)
            else #dosage
                genomat = convert_ds(Float64, vcffile; model = snpmodel,
                key="DS", impute=true, center=false, scale=false)
            end
            SnpArrays.makestream(pvalfile, "w") do io
                q = length(snpset)
                γ̂ = Vector{Float64}(undef, q)
                Z = zeros(size(fittednullmodel.mm, 1), q)
                if test == :score
                    ts = OrdinalMultinomialScoreTest(fittednullmodel.model, Z)
                    copyto!(ts.Z, @view(genomat[vcfrowinds, snpset]))
                    pval = polrtest(ts)
                    println(io, "The joint pvalue of snps indexed",
                        " at $(snpset) is $pval")
                elseif test == :lrt
                    nulldev = deviance(fittednullmodel.model)
                    Xaug = [fittednullmodel.model.X Z]
                    copyto!(@view(Xaug[:, fittednullmodel.model.p+1:end]), 
                            @view(genomat[vcfrowinds, snpset]))
                    altmodel = polr(Xaug, fittednullmodel.model.Y, 
                        fittednullmodel.model.link, solver, 
                        wts = fittednullmodel.model.wts)
                    copyto!(γ̂, 1, altmodel.β, fittednullmodel.model.p + 1, q)
                    l2normeffect = norm(γ̂)
                    pval = ccdf(Chisq(q), nulldev - deviance(altmodel))
                    println(io, "The l2norm of the effect size vector",
                    " is $l2normeffect and joint pvalue of snps indexed", 
                    " at $(snpset) is $pval")
                end
            end
        end
    else #analysistype == "gxe"
        isnothing(e) && 
            @error("GxE analysis indicated but not environmental variable keyword argument: `e` set.")
        # create holders for chromome, position, id 
        rec_chr = Array{Any, 1}(undef, 1)
        rec_pos = Array{Any, 1}(undef, 1)
        rec_ids = Array{Any, 1}(undef, 1)
        gholder = zeros(Union{Missing, Float64}, nsamples)
        Xaug = [fittednullmodel.model.X zeros(size(fittednullmodel.mm, 1))]
        Xaug2 = [fittednullmodel.model.X zeros(size(fittednullmodel.mm, 1), 2)] #or get Xaug to point to part of it

        # create array for environmental variable and testing 
        envvar = modelmatrix(FormulaTerm(fittednullmodel.mf.f.lhs, Term(Symbol(e))),
                 DataFrame(fittednullmodel.mf.data))
        testvec = Matrix{Float64}(undef, size(envvar))
        snpeffectnull = 0.0
        SnpArrays.makestream(pvalfile, "w") do io
            if test == :score 
                println(io, "chr,pos,snpid,snpeffectnull,pval")
            else 
                println(io, "chr,pos,snpid,snpeffectnull,snpeffectfull,GxEeffect,pval")
            end
            for j in eachindex(snpmask)
                if !snpmask[j] #skip snp, must read marker still. 
                    copy_gt!(gholder, reader; model = snpmodel, impute = true,
                    record_chr = rec_chr, record_pos = rec_pos, record_ids = rec_ids)
                    continue
                end
                if vcftype == :GT #genotype 
                    copy_gt!(gholder, reader; model = snpmodel, impute = true,
                    record_chr = rec_chr, record_pos = rec_pos, record_ids = rec_ids)
                else #dosage
                    copy_ds!(gholder, reader; model = snpmodel, impute = true,
                    record_chr = rec_chr, record_pos = rec_pos, record_ids = rec_ids)
                end
                copyto!(@view(Xaug[:, end]), @view(gholder[vcfrowinds]))
    
                if test == :score
                    copyto!(testvec, @view(Xaug[:, end]) .* envvar)
                    nm = polr(Xaug, fittednullmodel.model.Y, 
                        fittednullmodel.model.link, solver, wts = fittednullmodel.model.wts)
                    if Inf in nm.vcov #singular design matrix
                        snpeffectnull = 0.0
                        pval = 1.0
                    else
                        snpeffectnull = nm.β[end]
                        ts = OrdinalMultinomialScoreTest(nm, testvec)
                        pval = polrtest(ts)
                    end
                    println(io, "$(rec_chr[1]),$(rec_pos[1]),$(rec_ids[1][1]),$snpeffectnull,$pval")
                elseif test == :lrt
                    γ̂ = 0.0 # effect size for columns being tested
                    copyto!(@view(Xaug2[:, end - 1]), @view(Xaug[:, end]))
                    copyto!(@view(Xaug2[:, end]), @view(Xaug[:, end]) .*
                        envvar)
                    nm = polr(Xaug, fittednullmodel.model.Y, 
                    fittednullmodel.model.link, solver, wts = fittednullmodel.model.wts)
                    if Inf in nm.vcov #singular design matrix
                        snpeffectnull = 0.0
                        snpeffectfull = 0.0
                        γ̂ = 0.0
                        pval = 1.0
                    else
                        snpeffectnull = nm.β[end]
                        nulldev = deviance(nm)
                        altmodel = polr(Xaug2, fittednullmodel.model.Y, 
                            fittednullmodel.model.link, solver, 
                            wts = fittednullmodel.model.wts)
                        γ̂ = altmodel.β[end]
                        snpeffectfull = altmodel.β[end-1]
                        pval = ccdf(Chisq(1), nulldev - deviance(altmodel))
                    end
                    println(io, "$(rec_chr[1]),$(rec_pos[1]),$(rec_ids[1][1]),",
                        "$snpeffectnull,$snpeffectfull,$γ̂,$pval")
                end
            end
        end
    end
    close(reader)
    return fittednullmodel
end

# For BGEN Analysis
function ordinalgwas(
    fittednullmodel::StatsModels.TableRegressionModel,
    bgenfile::Union{AbstractString, IOStream}, # full path and bgen file name
    nsamples::Integer;          # number of samples in bed file
    analysistype::AbstractString = "singlesnp",
    testformula::FormulaTerm = fittednullmodel.mf.f.lhs ~ Term(:snp),
    test::Symbol = :score,
    pvalfile::Union{AbstractString, IOStream} = "ordinalgwas.pval.txt", 
    snpmodel::Union{Val{1}, Val{2}, Val{3}} = ADDITIVE_MODEL,
    snpinds::Union{Nothing, AbstractVector{<:Integer}} = nothing,
    bgenrowinds::AbstractVector{<:Integer} = 1:nsamples, # row indices for VCF array
    solver = NLoptSolver(algorithm=:LD_SLSQP, maxeval=4000),
    verbose::Bool = false,
    snpset::Union{Nothing, Integer, AbstractString, #for snpset analysis
        AbstractVector{<:Integer}} = nothing,
    e::Union{Nothing, AbstractString, Symbol} = nothing # for GxE analysis
    )

    # open BGEN file and get number of SNPs in file
    bgendata = Bgen(bgenfile)
    nsnps = n_variants(bgendata)
    bgen_iterator = iterator(bgendata)

    dosageholder = Vector{Float32}(undef, n_samples(bgendata))
    decompressed_length, _ = BGEN.check_decompressed_length(
        bgendata.io, first(bgen_iterator), bgendata.header)
    decompressed = Vector{UInt8}(undef, decompressed_length)

    # create SNP mask vector
    if snpinds == nothing
        snpmask = trues(nsnps)
    elseif eltype(snpinds) == Bool
        snpmask = snpinds
    else
        snpmask = falses(nsnps)
        snpmask[snpinds] .= true
    end

    analysistype = lowercase(analysistype)
    analysistype in ["singlesnp", "snpset", "gxe"] || error("Analysis type $analysis invalid option. 
    Available options are 'singlesnp', 'snpset' and 'gxe'.")


    # determine analysis type
    if analysistype == "singlesnp"
        # extra columns in design matrix to be tested
        testdf = DataFrame(fittednullmodel.mf.data) # TODO: not type stable here
        testdf[!, :snp] = zeros(size(fittednullmodel.mm, 1))
        Z = similar(modelmatrix(testformula, testdf))

        # create holder for dosage/snps 
        snpholder = zeros(Union{Missing, Float64}, size(fittednullmodel.mm, 1))

        # carry out score or LRT test SNP by SNP
        snponly = testformula.rhs == Term(:snp)
        SnpArrays.makestream(pvalfile, "w") do io
            if test == :score 
                println(io, "chr,pos,snpid,varid,pval")
                ts = OrdinalMultinomialScoreTest(fittednullmodel.model, Z)
            else 
                nulldev = deviance(fittednullmodel.model)
                Xaug = [fittednullmodel.model.X Z]
                q = size(Z, 2)
                γ̂ = Vector{Float64}(undef, q) # effect size for columns being tested
                if snponly
                    println(io, "chr,pos,snpid,varid,effect,stder,pval")
                else
                    print(io, "chr,pos,snpid,varid,")
                    for j in 1:q
                        print(io, "effect$j,")
                    end
                    println(io, "pval")
                end        
            end
            for (j, variant) in enumerate(bgen_iterator)
                if !snpmask[j] #skip snp
                    continue
                end
                minor_allele_dosage!(bgendata, variant; 
                    T = Float64, mean_impute = true, data = dosageholder, 
                    decompressed = decompressed)
                @views copyto!(snpholder, dosageholder[bgenrowinds])
                if test == :score
                    if snponly
                        copyto!(ts.Z, snpholder)
                    else # snp + other terms
                        copyto!(testdf[!, :snp], snpholder)
                        ts.Z[:] = modelmatrix(testformula, testdf)
                    end
                    pval = polrtest(ts)
                    println(io, "$(variant.chrom),$(variant.pos),$(variant.rsid),",
                    "$(variant.varid),$pval")
                elseif test == :lrt 
                    if snponly
                        copyto!(@view(Xaug[:, fittednullmodel.model.p+1]), 
                        snpholder)
                    else # snp + other terms
                        copyto!(testdf[!, :snp], snpholder)
                        Xaug[:, fittednullmodel.model.p+1:end] = modelmatrix(testformula, testdf)
                    end
                    altmodel = polr(Xaug, fittednullmodel.model.Y, 
                        fittednullmodel.model.link, solver, 
                        wts = fittednullmodel.model.wts)
                    copyto!(γ̂, 1, altmodel.β, fittednullmodel.model.p + 1, q)
                    pval = ccdf(Chisq(q), nulldev - deviance(altmodel))
                    stderr = stderror(altmodel)[fittednullmodel.model.npar + 1]
                    if snponly
                        println(io, "$(variant.chrom),$(variant.pos),$(variant.rsid),",
                        "$(variant.varid),$(γ̂[1]),$stderr,$pval")
                    else
                        print(io, "$(variant.chrom),$(variant.pos),$(variant.rsid),",
                        "$(variant.varid),")
                        for j in 1:q
                            print(io, "$(γ̂[j]),")
                        end
                        println(io, pval)
                    end
                end
            end
        end
    elseif analysistype == "snpset"
        # max size of a snpset length
        maxsnpset = 1

        #determine snpset
        if isa(snpset, Nothing)
            setlength = 1
            maxsnpset = 1
        elseif isa(snpset, AbstractString)
            isfile(snpset) || throw(ArgumentError("snpset file not found, 
            to specify a window replace snpset string with a window size"))
            #first column SNPset ID, second column SNP ID
            snpsetFile = CSV.read(snpset, DataFrame, header = [:snpset_id, :snp_id], delim = " ")
            maxsnpset = combine(groupby(snpsetFile, :snpset_id), :snp_id => length => :snpset_length) |> 
                x -> maximum(x.snpset_length)
            snpset_ids = unique(snpsetFile[!, :snpset_id])
            nSets = length(snpset_ids)
            setlength = 0
        elseif isa(snpset, Integer)
            setlength = snpset
            maxsnpset = snpset 
        else #abstract vector (boolean of true at indicies or range or indicies)
            setlength = -1
            maxsnpset = count(snpset .!= 0)
        end

        # create holders for chromome, position, id, dosage/gt
        snpholder = zeros(Union{Missing, Float64}, 
            size(fittednullmodel.mm, 1), maxsnpset)
        dosageholder = Vector{Float32}(undef, n_samples(bgendata))
        decompressed = Vector{UInt8}(undef, 3 * n_samples(bgendata) + 10)

        if setlength > 0 #single snp analysis or window
            Z = zeros(size(fittednullmodel.mm, 1), setlength) #
            q = setlength
            SnpArrays.makestream(pvalfile, "w") do io
                if test == :score
                    println(io, "startchr,startpos,startsnpid,startvarid,",
                        "endchr,endpos,endsnpid,endvarid,pval")
                    ts = OrdinalMultinomialScoreTest(fittednullmodel.model, Z)
                elseif test == :lrt 
                    println(io, "startchr,startpos,startsnpid,startvarid,",
                    "endchr,endpos,endsnpid,endvarid,l2normeffect,pval")
                    nulldev = deviance(fittednullmodel.model)
                    Xaug = [fittednullmodel.model.X Z]
                    γ̂ = Vector{Float64}(undef, setlength) # effect size for columns being tested
                end
                chrstart, posstart, rsidstart, varidstart = "", "", "", ""
                chrend, posend, rsidend, varidend = "", "", "", ""
                for j in 1:q:nsnps
                    endj = j + q - 1    
                    if endj >= nsnps
                        q = nsnps - j + 1
                        #length of Z will be different
                        snpholder = zeros(size(fittednullmodel.mm, 1), q)
                        if test == :score
                            Z = zeros(size(fittednullmodel.mm, 1), q)
                            ts = OrdinalMultinomialScoreTest(fittednullmodel.model, Z)
                        elseif test == :lrt 
                            Xaug = [fittednullmodel.model.X zeros(size(
                            fittednullmodel.mm, 1), q)]
                        end
                    end
                    for i in 1:q
                        variant = variant_by_index(bgendata, j + i - 1)
                        minor_allele_dosage!(bgendata, variant; 
                        T = Float64, mean_impute = true, data = dosageholder, 
                        decompressed = decompressed)
                        @views copyto!(snpholder[:, i], dosageholder[bgenrowinds])
                        if i == 1
                            chrstart = variant.chrom
                            posstart = variant.pos
                            rsidstart = variant.rsid
                            varidstart = variant.rsid
                        end
                        if i == q
                            chrend = variant.chrom
                            posend = variant.pos
                            rsidend = variant.rsid
                            varidend = variant.rsid
                        end
                    end
                    if test == :score
                        copyto!(ts.Z, snpholder)
                        pval = polrtest(ts)
                        println(io, "$chrstart,$posstart,$rsidstart,$varidstart,",
                        "$chrend,$posend,$rsidend,$varidend,$pval")
                    elseif test == :lrt 
                        copyto!(@view(Xaug[:, (fittednullmodel.model.p+1):end]), 
                        snpholder)
                        altmodel = polr(Xaug, fittednullmodel.model.Y, 
                        fittednullmodel.model.link, solver, 
                        wts = fittednullmodel.model.wts)
                        copyto!(γ̂, @view(altmodel.β[(fittednullmodel.model.p+1):end]))#, fittednullmodel.model.p + 1, setlength)
                        l2normeffect = norm(γ̂)
                        pval = ccdf(Chisq(q), nulldev - deviance(altmodel))
                        println(io, "$chrstart,$posstart,$rsidstart,$varidstart,",
                        "$chrend,$posend,$rsidend,$varidend,$l2normeffect,$pval")
                    end
                end
            end
        elseif setlength == 0 #snpset is defined by snpset file
            SnpArrays.makestream(pvalfile, "w") do io
                test == :score ? println(io, "snpsetid,nsnps,pval") : println(io, 
                    "snpsetid,nsnps,l2normeffect,pval")
                for j in eachindex(snpset_ids)
                    snpset_id = snpset_ids[j]
                    snpinds = findall(snpsetFile[!, :snpset_id] .== snpset_id)
                    q = length(snpinds)
                    Z = zeros(size(fittednullmodel.mm, 1), q)
                    for i in 1:q
                        variant = variant_by_index(bgendata, snpinds[i])
                        minor_allele_dosage!(bgendata, variant; 
                            T = Float64, mean_impute = true, data = dosageholder, 
                            decompressed = decompressed)
                        @views copyto!(Z[:, i], dosageholder[bgenrowinds])
                    end
                    if test == :score
                        ts = OrdinalMultinomialScoreTest(fittednullmodel.model, Z)
                        pval = polrtest(ts)
                        println(io, "$(snpset_id),$q,$pval")
                    elseif test == :lrt
                        γ̂ = Vector{Float64}(undef, q)
                        Xaug = [fittednullmodel.model.X Z]
                        nulldev = deviance(fittednullmodel.model)
                        altmodel = polr(Xaug, fittednullmodel.model.Y, 
                            fittednullmodel.model.link, solver, 
                            wts = fittednullmodel.model.wts)
                        copyto!(γ̂, 1, altmodel.β, fittednullmodel.model.p + 1, q)
                        l2normeffect = norm(γ̂)
                        pval = ccdf(Chisq(q), nulldev - deviance(altmodel))
                        println(io, "$(snpset_id),$q,$l2normeffect,$pval")
                    end
                end
            end
        else #setlength == -1 (testing just one set with specified snps in snpset)
            snpset = eltype(snpset) == Bool ? findall(snpset) : snpset 
            q = length(snpset)
            γ̂ = Vector{Float64}(undef, q)
            Z = zeros(size(fittednullmodel.mm, 1), q)
            for i in 1:length(snpset)
                variant = variant_by_index(bgendata, snpset[i])
                minor_allele_dosage!(bgendata, variant; 
                    T = Float64, mean_impute = true, data = dosageholder, 
                    decompressed = decompressed)
                @views copyto!(Z[:, i], dosageholder[bgenrowinds])
            end
            SnpArrays.makestream(pvalfile, "w") do io
                if test == :score
                    ts = OrdinalMultinomialScoreTest(fittednullmodel.model, Z)
                    pval = polrtest(ts)
                    println(io, "The joint pvalue of snps indexed",
                        " at $(snpset) is $pval")
                elseif test == :lrt
                    nulldev = deviance(fittednullmodel.model)
                    Xaug = [fittednullmodel.model.X Z]
                    altmodel = polr(Xaug, fittednullmodel.model.Y, 
                        fittednullmodel.model.link, solver, 
                        wts = fittednullmodel.model.wts)
                    copyto!(γ̂, 1, altmodel.β, fittednullmodel.model.p + 1, q)
                    l2normeffect = norm(γ̂)
                    pval = ccdf(Chisq(q), nulldev - deviance(altmodel))
                    println(io, "The l2norm of the effect size vector",
                    " is $l2normeffect and joint pvalue of snps indexed", 
                    " at $(snpset) is $pval")
                end
            end
        end
    else #analysistype == "gxe"
        isnothing(e) && 
            @error("GxE analysis indicated but not environmental variable keyword argument: `e` set.")

        Xaug = [fittednullmodel.model.X zeros(size(fittednullmodel.mm, 1))]
        Xaug2 = [fittednullmodel.model.X zeros(size(fittednullmodel.mm, 1), 2)] #or get Xaug to point to part of it

        # create array for environmental variable and testing 
        envvar = modelmatrix(FormulaTerm(fittednullmodel.mf.f.lhs, Term(Symbol(e))),
                 DataFrame(fittednullmodel.mf.data))
        testvec = Matrix{Float64}(undef, size(envvar))
        snpeffectnull = 0.0
        SnpArrays.makestream(pvalfile, "w") do io
            if test == :score 
                println(io, "chr,pos,snpid,varid,snpeffectnull,pval")
            else 
                println(io, "chr,pos,snpid,varid,snpeffectnull,snpeffectfull,GxEeffect,pval")
            end
            for (j, variant) in enumerate(bgen_iterator)
                if !snpmask[j] #skip snp, must read marker still. 
                    continue
                end
                minor_allele_dosage!(bgendata, variant; 
                    T = Float64, mean_impute = true, data = dosageholder, 
                    decompressed = decompressed)
                copyto!(@view(Xaug[:, end]), dosageholder[bgenrowinds])
    
                if test == :score
                    copyto!(testvec, @view(Xaug[:, end]) .* envvar)
                    nm = polr(Xaug, fittednullmodel.model.Y, 
                        fittednullmodel.model.link, solver, wts = fittednullmodel.model.wts)
                    if Inf in nm.vcov #singular design matrix
                        snpeffectnull = 0.0
                        pval = 1.0
                    else
                        snpeffectnull = nm.β[end]
                        ts = OrdinalMultinomialScoreTest(nm, testvec)
                        pval = polrtest(ts)
                    end
                    println(io, "$(variant.chrom),$(variant.pos),$(variant.rsid),",
                    "$(variant.varid),$snpeffectnull,$pval")
                elseif test == :lrt
                    γ̂ = 0.0 # effect size for columns being tested
                    copyto!(@view(Xaug2[:, end - 1]), @view(Xaug[:, end]))
                    copyto!(@view(Xaug2[:, end]), @view(Xaug[:, end]) .*
                        envvar)
                    nm = polr(Xaug, fittednullmodel.model.Y, 
                    fittednullmodel.model.link, solver, wts = fittednullmodel.model.wts)
                    if Inf in nm.vcov #singular design matrix
                        snpeffectnull = 0.0
                        snpeffectfull = 0.0
                        γ̂ = 0.0
                        pval = 1.0
                    else
                        snpeffectnull = nm.β[end]
                        nulldev = deviance(nm)
                        altmodel = polr(Xaug2, fittednullmodel.model.Y, 
                            fittednullmodel.model.link, solver, 
                            wts = fittednullmodel.model.wts)
                        γ̂ = altmodel.β[end]
                        snpeffectfull = altmodel.β[end-1]
                        pval = ccdf(Chisq(1), nulldev - deviance(altmodel))
                    end
                    println(io, "$(variant.chrom),$(variant.pos),$(variant.rsid),",
                        "$(variant.varid),$snpeffectnull,$snpeffectfull,$γ̂,$pval")
                end
            end
        end
    end
    return fittednullmodel
end

# VCFTools uses different coding for additive, dominant, recessive models than SnpArrays
modelingdict = Dict(
    Val{1}() => :additive,
    Val{2}() => :dominant,
    Val{3}() => :recessive
    )