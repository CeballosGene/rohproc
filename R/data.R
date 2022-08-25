#' ROH in the 1K Genomes.
#'
#' This dataset contains the PLINKs .hom file obtained using 1K Genomes array data.
#' PLINKs parameters were set to:
#' homozyg-window-snp 30.
#' homozyg-window-het 1.
#' homozyg-snp 30.
#' homozyg-kb 300.
#' homozyg-density 30.
#' @format A data frame with 754450 rows and 13 variables:
#' \describe{
#'   \item{IID}{individual id}
#'   \item{FID}{family id}
#'   \item{PHE}{Phenotype not present}
#'   \item{CHR}{Chromosome number}
#'   \item{SNP1}{rs of the SNP in starting possition}
#'   \item{SNP2}{rs of the SNP in ending possition}
#'   \item{POS1}{starting of the ROH}
#'   \item{POS2}{ending of the ROH}
#'   \item{KB}{ROH length in KB}
#'   \item{NSNP}{Number of SNP present in the ROH}
#'   \item{DENSITY}{Densito of SNP}
#'   \item{PHOM}{}
#'   \item{PHET}{}
#' }
#' @source \url{http://https://www.internationalgenome.org}
"KGenomes_hom"

#' Populations 1K Genomes.
#'
#' This dataset contains the classification in populations for all the individuals in the 1K Genomes dataset:
#' @format A data frame with 19 rows and 2 variables:
#' \describe{
#'   \item{IID}{individual id}
#'   \item{pop{population}
#' }
#' @source \url{http://https://www.internationalgenome.org}
"KGenomes_pops"

#' Continents in 1K Genomes.
#'
#' This dataset contains the classification in continents for all the populations in the 1K Genomes dataset:
#' @format A data frame with 19 rows and 2 variables:
#' \describe{
#'   \item{pop}{population}
#'   \item{cont}{continent}
#' }
#' @source \url{http://https://www.internationalgenome.org}
"KGenomes_cont"
