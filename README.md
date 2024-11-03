### Pore size eDNA

**About**

This repository contains the environments, pipelines, scripts, notebooks and auxiliary files belonging to the publication about eDNA pore size on seawater danish samples.

*Filter pore size influences taxonomic composition of retained eDNA from seawater samples – evidence from shotgun sequencing*

**Usage**

The raw data can be analyzed, on a HPC cluster using Slurm, by:

    - Installing the available environments (folder: envs) on conda/mamba
    - Installing ROBITools and ROBITaxonomy inside the ROBITools environment (details below)
    - Run the pipelines:
        - eDNA.first.py: Data merging and splitting followed by Blastn mapping and filter.
        - eDNA.second.py: LCA of filtered results, merging of results into unique tables per sample and downstream analysis.
    - Once the data is ready, the notebooks could be run to summarize results and produce the main and supplementary figures.
    - Environment usage:
        - eDNA: eDNA.first.py
        - ROBITools and TaxizeDB: eDNA.second.py
        - Jup-eDNA: jupyter notebooks

The pipeline was splitted in two to improve the tracking of scheduled jobs (There are many thousands of Blastn jobs).

**Installing ROBITools/Taxonomy**

```
# Install in ROBITools Environment
# devtools::install_gitlab("obitools/ROBITaxonomy", host = "https://git.metabarcoding.org/")
# install.packages("igraph")
# devtools::install_gitlab("obitools/ROBITools", host = "https://git.metabarcoding.org/")
```