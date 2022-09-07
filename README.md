rohproc R Package
================
Francisco C. Ceballos
2022-09-07

-   <a href="#package-overview" id="toc-package-overview">Package
    Overview</a>
    -   <a href="#data-avaliable-with-this-package"
        id="toc-data-avaliable-with-this-package">Data avaliable with this
        package</a>
-   <a href="#processing-hom-file" id="toc-processing-hom-file">PROCESSING
    .HOM FILE</a>
-   <a href="#population-genetics" id="toc-population-genetics">POPULATION
    GENETICS</a>
-   <a href="#roh-islands-and-regions-of-heterozygosity"
    id="toc-roh-islands-and-regions-of-heterozygosity">ROH islands and
    Regions of Heterozygosity.</a>

## Package Overview

The analysis of autozygosity through the Runs of Homozygosity (ROH)
–contiguous regions of the genome where an individual is homozygous
across all sites– is an alternative and relatively new approach to
answer many important questions in modern genetics, including population
history and the genetic architecture of complex traits.

This package provides with several useful functions to represent and
analyze PLINK’s ROH outcomes. The following script in bash is
recommended to search of ROH in humans. The following PLINK’s conditions
are optimal for ROH searching in humans when a minimun of 1.5M SNP is
used.

<pre><code>plink  \
  --bfile path/to/bfile \
  --homozyg-window-snp 30 \
  --homozyg-window-het 1 \
  --homozyg-snp 30 \
  --homozyg-kb 300 \
  --homozyg-density 30 \
  --out path/to/folder
</code></pre>

This bash script will generate a .hom file with all the information
needed. This package takes that tabular file and process it to obtain
different summary outcomes, figures and other features.

### Data avaliable with this package

This package includes .hom file with the PLINK’s ROH analysis of 19
populations belonging to the [1000 Genomes
Project](https://www.internationalgenome.org/).

## PROCESSING .HOM FILE

First three functions avaliable in this package are useful to first
process the PLINK’s .hom file.

## POPULATION GENETICS

Two important data representations for population genetics insights:
\> - *Number of ROH vs Sum of ROH*

> -   *FIS vs FROH*

## ROH islands and Regions of Heterozygosity.
