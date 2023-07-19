# Isosceles 

Isosceles (**Iso**forms from **S**ingle-**Ce**ll; **L**ong-read **E**xpression 
**S**uite) is an R package dedicated to transcript detection and quantification 
from ONT reads, supporting both bulk RNA-Seq and scRNA-Seq technologies.

Isosceles can be installed using the following commands:
```r
# Install Bioconductor dependencies (basic installation)
BiocManager::install(c("scuttle", "GenomicFeatures", "BSgenome"))
# Install Bioconductor dependencies (building the vignette)
BiocManager::install(c("dittoSeq", "ggbio", "BiocStyle"))
# Install Isosceles
devtools::install_github('timbitz/Isosceles')
```

Load the Isosceles package:
```r
library(Isosceles)
```

You can follow along with our vignette [here](https://timbitz.github.io/Isosceles/inst/vignette/Isosceles.html)!
