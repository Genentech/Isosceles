# Isosceles 

Isosceles (**Iso**forms from **S**ingle-**Ce**ll; **L**ong-read **E**xpression 
**S**uite) is an R package dedicated to transcript detection and quantification 
from long reads, supporting both bulk RNA-Seq and scRNA-Seq technologies.

<p align="center">
  <img src="docs/Isosceles_header.gif" width="600">
</p>


Kabza M., Ritter A., Byrne A., Sereti K., Le D., Stephenson W., Sterne-Weiler T. Accurate long-read transcript discovery and quantification at single-cell, pseudo-bulk and bulk resolution with Isosceles. _Nat Commun_ **15**, 7316 (2024). https://doi.org/10.1038/s41467-024-51584-3

## Installation

Isosceles can be installed using the following commands:
```r
install.packages(c("BiocManager", "devtools"))
BiocManager::install(c("scran", "scater", "uwot", "dittoSeq", "DEXSeq", 
                       "Nebulosa", "ggbio", "BiocStyle"))
devtools::install_github("Genentech/Isosceles", dependencies = TRUE, upgrade = TRUE,
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

You can follow along with our vignettes ([Introduction to the Isosceles package](https://genentech.github.io/Isosceles/docs/Isosceles.html), [Mouse E18 brain data analysis](https://genentech.github.io/Isosceles/docs/Mouse_E18_brain_analysis.html))
or the [reference manual](https://github.com/Genentech/Isosceles/blob/devel/docs/Isosceles.pdf)!

## Best practices

  * We recommend [minimap2](https://github.com/lh3/minimap2) for all long-read alignments.
  * **Isosceles doesn't perform post-hoc splice junction correction, so is critical to run minimap2 with the '\-\-junc-bed' flag.** The intron position BED file required by it can be easily created using the `gtf_to_intron_bed` function.
  * The default settings of de novo transcript detection used by the `bam_to_tcc` function should work well for most expected read depths across eukaryotic transcriptomes, but for the analysis of spike-in data, such as SIRVs, we recommend increasing the read count threshold (the *min_read_count* argument) to a higher value (e.g. 50).

## Troubleshooting

You can check if the package works correctly by running its unit tests:
```r
testthat::test_package("Isosceles")
```

In case of any problems, we recommend using the Isosceles Singularity image you
can download from [Zenodo](https://zenodo.org/doi/10.5281/zenodo.8180648)
or installing the package in a
[Docker container using a Bioconductor image](https://www.bioconductor.org/help/docker).
