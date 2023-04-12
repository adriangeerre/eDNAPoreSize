## Pore Size eDNA pipeline

#### Installation

```
mamba create -n eDNA_mamba -c bioconda -c conda-forge -c r -c gwforg bbmap blast fastqc seqtk gwf=1.7.2 r-vegan r-tidyverse r-plyr r-reshape2 r-ggalt r-ggpubr r-ggrepel r-taxizedb bioconductor-complexheatmap bioconductor-edger
```


```
conda activate eDNA_mamba
R
    install.packages("optparse")
```