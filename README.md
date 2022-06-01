[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
![snakemaker](https://github.com/AnimalGenomicsETH/bovine-assembly/workflows/snakemaker/badge.svg?branch=master)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.6503779.svg)](https://doi.org/10.5281/zenodo.6503779)
# Bovine pangenomes


## Genome assembly

This repository contains pipelines and scripts used to assemble, validate, and compare bovine genomes using ONT and HiFi data for three trios
  - Original Braunvieh x Original Braunvieh
  - Nellore x Brown Swiss
  - gaur x Piedmontese

For multiple current long read assemblers assemblers like hifiasm, canu, Shasta, Flye, etc.


## Pangenome analysis

There are also multiple pipelines and scripts used to construct different pangenomes with minigraph and subsequent analysis of the SVs like repeatitive content and accuracy scores.

## Citation

Details on the results can be found in our [publication](https://www.nature.com/articles/s41467-022-30680-2), or the assemblies hosted [here](https://doi.org/10.5281/ZENODO.5906579). 

> **Structural variant-based pangenome construction has low sensitivity to variability of haplotype-resolved bovine assemblies**
> 
> Alexander S. Leonard, Danang Crysnanto, Zih-Hua Fang, Michael P Heaton, Brian L. Vander Ley, Carolina Herrera, Heinrich Bollwein, Derek M. Bickhart, Kristen L. Kuhn, Timothy PL. Smith, Benjamin D. Rosen, Hubert Pausch
> 
> _Nat Commun_ 13, 3012 (2022). https://doi.org/10.1038/s41467-022-30680-2
> 


## Usage
Many of the parameters are tuned to run for our data and on the ETH Euler cluster, using for example a forked version of the LSF snakemake [profile](https://github.com/AnimalGenomicsETH/snakemake_lsf), so it may take some modifying to work smoothly in different contexts. Many tools are assumed to be available in $PATH, but all are freely available.
