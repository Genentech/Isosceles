# Isosceles 

Isosceles (**Iso**forms from **S**ingle-**Ce**ll; **L**ong-read **E**xpression 
**S**uite) is an R package dedicated to transcript detection and quantification 
from long reads, supporting both bulk RNA-Seq and scRNA-Seq technologies.

<p align="center">
  <img src="docs/Isosceles_header.gif" width="600">
</p>

Isosceles can be installed using the following commands:
```r
BiocManager::install(c("scran", "scater", "uwot", "dittoSeq", "DEXSeq", 
                       "Nebulosa", "ggbio", "BiocStyle"))
devtools::install_github("timbitz/Isosceles", dependencies = TRUE, upgrade = TRUE)
```

Load the Isosceles package:
```r
library(Isosceles)
```

You can follow along with our vignettes ([Introduction to the Isosceles package](https://timbitz.github.io/Isosceles/docs/Isosceles.html))
or the [reference manual](https://github.com/timbitz/Isosceles/blob/devel/docs/Isosceles.pdf)!
