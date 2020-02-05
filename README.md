# OrdinalGWAS

| **Documentation** | **Build Status** | **Code Coverage**  |
|-------------------|------------------|--------------------|
| [![](https://img.shields.io/badge/docs-stable-blue.svg)](https://OpenMendel.github.io/OrdinalGWAS.jl/stable) [![](https://img.shields.io/badge/docs-latest-blue.svg)](https://OpenMendel.github.io/OrdinalGWAS.jl/latest) | [![Build Status](https://travis-ci.org/OpenMendel/OrdinalGWAS.jl.svg?branch=master)](https://travis-ci.org/OpenMendel/OrdinalGWAS.jl) [![Build status](https://ci.appveyor.com/api/projects/status/9qkvj4pbw8cwxa9l/branch/master?svg=true)](https://ci.appveyor.com/project/Hua-Zhou/OrdinalGWAS-jl/branch/master) | [![Coverage Status](https://coveralls.io/repos/github/OpenMendel/OrdinalGWAS.jl/badge.svg?branch=master)](https://coveralls.io/github/OpenMendel/OrdinalGWAS.jl?branch=master) [![codecov](https://codecov.io/gh/OpenMendel/PolrGWAS.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/OpenMendel/OrdinalGWAS.jl) |  


OrdinalGWAS.jl is a Julia package for performing genome-wide association studies (GWAS) for ordered categorical phenotypes. It is useful when the phenotype takes ordered discrete values, e.g., disease status (undiagnosed, pre-disease, mild, moderate, severe). It runs efficiently and scales well to very large datasets. 

OrdinalGWAS.jl supports Julia v1.0 or later. See documentation for usage.  
[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://OpenMendel.github.io/OrdinalGWAS.jl/stable) [![](https://img.shields.io/badge/docs-latest-blue.svg)](https://OpenMendel.github.io/OrdinalGWAS.jl/latest)

OrdinalGWAS.jl is not yet registered. It requires [SnpArrays.jl](https://github.com/OpenMendel/SnpArrays.jl) and [OrdinalMultinomialModels.jl](https://github.com/OpenMendel/OrdinalMultinomialModels.jl) which are also not yet registered, so it will require the following steps to install. 

```{julia}
pkg> add https://github.com/OpenMendel/SnpArrays.jl.git

pkg> add https://github.com/OpenMendel/OrdinalMultinomialModels.jl.git

pkg> add https://github.com/OpenMendel/OrdinalGWAS.jl.git
```

## Citation

The methods and applications of this software package are detailed in the following publication:

*German CA, Sinsheimer JS, Klimentidis YC, Zhou H, Zhou JJ. Ordered multinomial regression for genetic association analysis of ordinal phenotypes at Biobank scale. Genetic Epidemiology. 2019; 1-13. https://doi.org/10.1002/gepi.22276*

If you use [OpenMendel](https://openmendel.github.io) analysis packages in your research, please cite the following reference in the resulting publications:

*OPENMENDEL: a cooperative programming project for statistical genetics. Zhou H, Sinsheimer JS, Bates DM, Chu BB, German CA, Ji SS, Keys KL, Kim J, Ko S, Mosher GD, Papp JC, Sobel EM, Zhai J, Zhou JJ, Lange K. Hum Genet. 2019 Mar 26. doi: 10.1007/s00439-019-02001-z. [Epub ahead of print] PMID: 30915546*

## Acknowledgments

This project is supported by the National Institutes of Health under NIGMS awards R01GM053275 and R25GM103774 and NHGRI award R01HG006139.



