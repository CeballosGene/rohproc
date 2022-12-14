rohproc R Package
================
Francisco C. Ceballos
2022-09-19

-   <a href="#package-overview" id="toc-package-overview">Package
    Overview</a>
    -   <a href="#data-avaliable-with-this-package"
        id="toc-data-avaliable-with-this-package">Data avaliable with this
        package</a>
-   <a href="#processing-hom-files" id="toc-processing-hom-files">Processing
    .hom files</a>
-   <a href="#population-genetics" id="toc-population-genetics">Population
    Genetics.</a>
    -   <a href="#roh-distribution" id="toc-roh-distribution">ROH
        distribution</a>
    -   <a href="#searching-for-the-origin-of-the-inbreeding"
        id="toc-searching-for-the-origin-of-the-inbreeding">Searching for the
        Origin of the inbreeding</a>
-   <a href="#roh-islands-and-regions-of-heterozygosity"
    id="toc-roh-islands-and-regions-of-heterozygosity">ROH islands and
    Regions of Heterozygosity.</a>
    -   <a href="#roh-islands-rohi" id="toc-roh-islands-rohi">ROH islands
        (ROHi)</a>
    -   <a href="#regions-of-heterozygosity-rhz"
        id="toc-regions-of-heterozygosity-rhz">Regions of heterozygosity
        (RHZ)</a>

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

## Processing .hom files

First three functions avaliable in this package are useful to first
process the PLINK’s .hom file.

## Population Genetics.

### ROH distribution

### Searching for the Origin of the inbreeding

Two distinct and independent biological scenarios can increase
homozygosity in natural populations: cultural consanguinity and genetic
drift in isolated populations. These two different sources were defined
in classical population genetics as systematic inbreeding (denoted by
FIS) and panmictic nbreeding (denoted by FST), respectively. Total
inbreeding, denoted by FIT, is defined by (1-FIT) = (1-FIS) (1-FST)
Panmictic inbreeding occurs in isolated populations, when individuals
randomly mate within their own group, with no immigration. Population
isolation can be cultural, a consequence of geographical barriers or
because of sedentary behavior. Importantly, isolation by itself does not
create genomic autozygosity, except when the effective population size
(Ne) is small and genetic drift has the strength to remove genetic
variability. On the other hand, systematic inbreeding, or cultural
consanguinity, has the effect of reducing heterozygosity relative to the
expectation under Hardy-Weinberg equilibrium independent of Ne, and thus
increasing FIS. High consanguinity (and consequent high FIS), and
genetic drift by isolation coupled with low Ne (and consequent high FST)
are two independent and non-mutually-exclusive phenomena that can
increase overall autozygosity (FIT) in a population. Here we also note
that the term ‘’endogamy’’ is generally used to describe population
isolation, although the term is sometimes used to refer to
consanguinity.

> -   *Number of ROH vs Sum of ROH* Wwe may assess the origins of the
>     autozygosity using comparisons of NROH\>1.5Mb and SROH\>1.5Mb.
>     Namely, if an individual displays excess of SROH\>1.5Mb relative
>     to NROH\>1.5Mb in comparison to non-consanguineous individuals,
>     this suggests autozygosity by consanguinity. If both SROH\>1.5Mb
>     and NROH\>1.5Mb are high, this suggests drift-driven autozygosity.
>     This approach does not provide a quantitative estimation of the
>     elative contributions of drift versus consanguinity, but only a
>     qualitative assessment. Nevertheless, it is a powerful approach:
>     its inferences on autozygosity patterns in modern human
>     populations’ genomes are consistent with known ethnographic data
>     about consanguineous traditions, and about population size and
>     isolation.Simulations of the number and sum of ROHs, for ROH \>
>     1.5 Mb, calculated for the offspring of different consanguineous
>     matings can be shown, along with the ancient and modern samples.
>     Points with different colors designate offspring of different
>     consanguineous mating: second cousin (green), first cousin
>     (yellow), avuncular (uncle-niece, aunt-nephew, double first
>     cousin) (orange), incest (brother-sister, parent-offspring) (red).
>     Five thousand simulations are represented for each consanguineous
>     mating. Note that this simulation does not include drift, but the
>     degree of right shift can be projected to cases where there exists
>     a non-zero level of autozygosity due to drift

![N_vs_Sum](C:/Users/Rembukai/Desktop/R_scripts/ROH_pack/imagenes/n_vs_s.jpg)

> -   *FIS vs FROH*

## ROH islands and Regions of Heterozygosity.

### ROH islands (ROHi)

ROH islands are defined as regions in the genome where the proportion of
individuals of a population deviates from the expected under a binomial
distribution. These regions have been found to be enriched with protein
coding genes under selection. To search for ROHi, a sliding window of
100 kb was used. In every 100 kb genomic window, the number of subjects
with ROH was obtained and a binomial test was applied (threshold for
significance established at p\<2x10-6, corresponding to an adjustment
for 25,000 windows).

### Regions of heterozygosity (RHZ)

RHZ are regions in the genome where no individual in the population has
a ROH. In order to only identify informative heterozygous haplotypes,
regions that have anomalous, unstructured, high signal/read counts in
next generation sequence experiments were removed. These 226 regions,
called ultra-high signal artefact regions, include high mapability
islands, low mapability islands, satellite repeats, centromere regions,
snRNA and telomeric regions (Consortium EP 2012).
