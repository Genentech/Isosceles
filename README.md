# Isosceles 

Isosceles (**Iso**forms from **S**ingle-**Ce**ll; **L**ong-read **E**xpression 
**S**uite) is an R package dedicated to transcript detection and quantification 
from long reads, supporting both bulk RNA-Seq and scRNA-Seq technologies.

<p align="center">
  <img src="docs/Isosceles_header.gif" width="600">
</p>

Preprint: https://www.biorxiv.org/content/10.1101/2023.11.30.566884

## Installation

Isosceles can be installed using the following commands:
```r
install.packages(c("BiocManager", "devtools"))
BiocManager::install(c("scran", "scater", "uwot", "dittoSeq", "DEXSeq", 
                       "Nebulosa", "ggbio", "BiocStyle"))
devtools::install_github("timbitz/Isosceles", dependencies = TRUE, upgrade = TRUE,
                         INSTALL_opts = "--install-tests")
```

We found that some versions of Isosceles' dependencies don't work together well, 
which might cause problems with testing the package or building the vignettes. 
If you encounter such issues, re-installing certain packages might be helpful: 
```r
install.packages("irlba") 
devtools::install_github("powellgenomicslab/Nebulosa", upgrade = FALSE) 
```

Load the Isosceles package:
```r
library(Isosceles)
```

## Usage

You can follow along with our vignettes ([Introduction to the Isosceles package](https://timbitz.github.io/Isosceles/docs/Isosceles.html), [Mouse E18 brain data analysis](https://timbitz.github.io/Isosceles/docs/Mouse_E18_brain_analysis.html))
or the [reference manual](https://github.com/timbitz/Isosceles/blob/devel/docs/Isosceles.pdf)!

## Troubleshooting

You can check if the package works correctly by running its unit tests:
```r
testthat::test_package("Isosceles")
```

In case of any problems, we recommend using the Isosceles Singularity image you
can download from [Zenodo](https://zenodo.org/doi/10.5281/zenodo.8180648)
or installing the package in a
[Docker container using a Bioconductor image](https://www.bioconductor.org/help/docker).
