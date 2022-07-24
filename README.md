# OrdinalGWAS

| **Documentation** | **Build Status** | **Code Coverage**  |
|-------------------|------------------|--------------------|
| [![](https://img.shields.io/badge/docs-stable-blue.svg)](https://OpenMendel.github.io/OrdinalGWAS.jl/stable) [![](https://img.shields.io/badge/docs-latest-blue.svg)](https://OpenMendel.github.io/OrdinalGWAS.jl/latest) | [![build Actions Status](https://github.com/OpenMendel/OrdinalGWAS.jl/workflows/CI/badge.svg)](https://github.com/OpenMendel/OrdinalGWAS.jl/actions) | [![Coverage Status](https://coveralls.io/repos/github/OpenMendel/OrdinalGWAS.jl/badge.svg?branch=master)](https://coveralls.io/github/OpenMendel/OrdinalGWAS.jl?branch=master) [![codecov](https://codecov.io/gh/OpenMendel/OrdinalGWAS.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/OpenMendel/OrdinalGWAS.jl) |  


OrdinalGWAS.jl is a Julia package for performing genome-wide association studies (GWAS) for ordered categorical phenotypes. It is useful when the phenotype takes ordered discrete values, e.g., disease status (undiagnosed, pre-disease, mild, moderate, severe). It runs efficiently and scales well to very large datasets. The package currently supports [PLINK](https://zzz.bwh.harvard.edu/plink/), [VCF](https://en.wikipedia.org/wiki/Variant_Call_Format) (both dosage and genotype data) file formats, and [BGEN](https://www.well.ox.ac.uk/~gav/bgen_format/) file formats. We plan to add [PGEN](https://www.cog-genomics.org/plink/2.0/formats#pgen) support in the future. 

OrdinalGWAS.jl supports Julia v1.6 or later. See the [documentation](https://openmendel.github.io/OrdinalGWAS.jl/latest/) for usage.  
[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://OpenMendel.github.io/OrdinalGWAS.jl/stable) [![](https://img.shields.io/badge/docs-latest-blue.svg)](https://OpenMendel.github.io/OrdinalGWAS.jl/latest)

OrdinalGWAS.jl needs the following steps to install. 

```{julia}
using Pkg
pkg"add OrdinalGWAS"
```

## Citation

The methods and applications of this software package are detailed in the following publication:

*German CA, Sinsheimer JS, Klimentidis YC, Zhou H, Zhou JJ. Ordered multinomial regression for genetic association analysis of ordinal phenotypes at Biobank scale. Genet Epidemiol. 2020 Apr;44(3):248-260. doi: 10.1002/gepi.22276. Epub 2019 Dec 26. PMID: 31879980; PMCID: [PMC8256450](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8256450/).*

If you use [OpenMendel](https://openmendel.github.io) analysis packages in your research, please cite the following reference in the resulting publications:

*Zhou H, Sinsheimer JS, Bates DM, Chu BB, German CA, Ji SS, Keys KL, Kim J, Ko S, Mosher GD, Papp JC, Sobel EM, Zhai J, Zhou JJ, Lange K. OPENMENDEL: a cooperative programming project for statistical genetics. Hum Genet. 2020 Jan;139(1):61-71. doi: 10.1007/s00439-019-02001-z. Epub 2019 Mar 26. PMID: 30915546; PMCID: [PMC6763373](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6763373/).*

## Acknowledgments

This project has been supported by the National Institutes of Health under awards R01GM053275, R01HG006139, R25GM103774, and 1R25HG011845.
