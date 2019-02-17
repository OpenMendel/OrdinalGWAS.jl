#!/usr/local/bin/julia
#
# This script demonstrates how to submit multiple OrdinalGWAS runs from multiple sets of
# Plink files on UCLA Hoffman2 cluster. It assumes that a demo data is available by
# running `julia cluster_preparedata.jl` at current folder.
#

using OrdinalGWAS, Serialization

# Step 1: fit null model and save result to file `fittednullmodel.jls`
nm = ordinalgwas(@formula(trait ~ sex), "covariate.txt", nothing)
open("fittednullmodel.jls", "w") do io
    Serialization.serialize(io, nm)
end

# Step 2: GWAS for each chromosome
for chr in 1:23
    println("submit job for chromosome=$chr")
    jcode = "using OrdinalGWAS, Serialization;
    nm = open(deserialize, \"fittednullmodel.jls\");
    bedfile = \"hapmap3.chr.\" * string($chr) * \".bed.gz\";
    bimfile = \"hapmap3.chr.\" * string($chr) * \".bim.gz\";
    pvalfile = \"hapmap3.chr.\" * string($chr) * \".pval.txt\";
    ordinalgwas(nm, bedfile, bimfile, 324; pvalfile=pvalfile);"
    # prepare sh file for qsub
    open("tmp.sh", "w") do io
        println(io, "#!/bin/bash")
        println(io, "#\$ -cwd")
        println(io, "# error = Merged with joblog")
        println(io, "#\$ -o joblog.\$JOB_ID")
        println(io, "#\$ -j y")
        println(io, "#\$ -l h_rt=0:30:00,h_data=2G") # request runtime and memory
        println(io, "#\$ -pe shared 2") # request # shared-memory nodes
        println(io, "# Email address to notify")
        println(io, "#\$ -M \$USER@mail")
        println(io, "# Notify when")
        println(io, "#\$ -m a")
        println(io)
        println(io, "# load the job environment:")
        println(io, ". /u/local/Modules/default/init/modules.sh")
        println(io, "module load julia/1.0.1") # available Julia version
        println(io)
        println(io, "# run julia code")
        println(io, "julia -e '$jcode' > output.\$JOB_ID 2>&1")
    end
    # submit job
    run(`qsub tmp.sh`)
end
