[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
![snakemaker](https://github.com/AnimalGenomicsETH/bovine-assembly/workflows/snakemaker/badge.svg?branch=master)

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

Details on the results can be found in our [manuscript](https://www.biorxiv.org/content/10.1101/2021.11.02.466900v2), or the assemblies hosted [here](https://doi.org/10.5281/ZENODO.5906579). 

> **Structural variant-based pangenome construction has low sensitivity to variability of haplotype-resolved bovine assemblies**
> 
> Alexander S. Leonard, Danang Crysnanto, Zih-Hua Fang, Michael P Heaton, Brian L. Vander Ley, Carolina Herrera, Heinrich Bollwein, Derek M. Bickhart, Kristen L. Kuhn, Timothy PL. Smith, Benjamin D. Rosen, Hubert Pausch
> 
> bioRxiv 2021.11.02.466900; doi: https://doi.org/10.1101/2021.11.02.466900


## Usage
Many of the parameters are tuned to run for our data and on the ETH Euler cluster, using for example a forked version of the LSF snakemake [profile](https://github.com/AnimalGenomicsETH/snakemake_lsf), so it may take some modifying to work smoothly in different contexts. Many tools are assumed to be available in $PATH, but all are freely available.
