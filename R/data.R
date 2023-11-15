#' ROH in the Human Genome Diversity Panel.
#'
#' This dataset contains the PLINKs .hom file obtained using Human Genome Diversity Panel.
#' PLINKs parameters were set to:
#' homozyg-window-snp 30.
#' homozyg-window-het 1.
#' homozyg-snp 30.
#' homozyg-kb 300.
#' homozyg-density 30.
#' @format ## HGDP_hom
#' A data frame with 612139 rows and 13 variables:
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
#' @source \url{https://cephb.fr/en/hgdp_panel.php}
"HGDP_hom"

#' Populations Human Genome Diversity Panel.
#'
#' This dataset contains the classification in populations for all the individuals in the HGD dataset:
#' @format A data frame with 57 rows and 2 variables:
#' \describe{
#'   \item{IID}{individual id}
#'   \item{pop}{population}
#' }
#' @source \url{https://cephb.fr/en/hgdp_panel.php}
"HGDP_pops"

#' Continents in Human Genome Diversity Panel.
#'
#' This dataset contains the classification in continents for all the populations in the HGD dataset:
#' @format A data frame with 57 rows and 2 variables:
#' \describe{
#'   \item{pop}{population}
#'   \item{cont}{continent}
#' }
#' @source \url{https://cephb.fr/en/hgdp_panel.php}
"HGDP_cont"

#' FIS in Human Genome Diversity Panel.
#'
#' This dataset contains the average FIS estimates for each individual in the HGD dataset:
#' @format A data frame with 936 rows and 2 variables:
#' \describe{
#'   \item{IID}{individual id}
#'   \item{Fis}{FIS estimate obtained using PLINK's --het flag}
#' }
#' @source \url{https://cephb.fr/en/hgdp_panel.php}
"HGDP_het"

#' Cluster Africa Human Genome Diversity Panel.
#'
#' This dataset contains the PLINKs .hom file obtained using Human Genome Diversity Panel. Cluster Africa.
#' PLINKs parameters were set to:
#' homozyg-window-snp 30.
#' homozyg-window-het 1.
#' homozyg-snp 30.
#' homozyg-kb 300.
#' homozyg-density 30.
#' @format A data frame with 26017 rows and 16 variables:
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
#'   \item{POP}{}
#'   \item{CLUSTER}{}
#'   \item{NAME_CLUS}{}
#' }
#' @source \url{https://cephb.fr/en/hgdp_panel.php}
"HGDP_cl_Africa"

#' Cluster East Asia Human Genome Diversity Panel.
#'
#' This dataset contains the PLINKs .hom file obtained using Human Genome Diversity Panel. Cluster Africa.
#' PLINKs parameters were set to:
#' homozyg-window-snp 30.
#' homozyg-window-het 1.
#' homozyg-snp 30.
#' homozyg-kb 300.
#' homozyg-density 30.
#' @format A data frame with 147196 rows and 16 variables:
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
#'   \item{POP}{}
#'   \item{CLUSTER}{}
#'   \item{NAME_CLUS}{}
#' }
#' @source \url{https://cephb.fr/en/hgdp_panel.php}
"HGDP_cl_EastAsia"

#' Cluster Europe Human Genome Diversity Panel.
#'
#' This dataset contains the PLINKs .hom file obtained using Human Genome Diversity Panel. Cluster Africa.
#' PLINKs parameters were set to:
#' homozyg-window-snp 30.
#' homozyg-window-het 1.
#' homozyg-snp 30.
#' homozyg-kb 300.
#' homozyg-density 30.
#' @format A data frame with 99960 rows and 16 variables:
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
#'   \item{POP}{}
#'   \item{CLUSTER}{}
#'   \item{NAME_CLUS}{}
#' }
#' @source \url{https://cephb.fr/en/hgdp_panel.php}
"HGDP_cl_Europe"

#' Cluster Middle East Human Genome Diversity Panel.
#'
#' This dataset contains the PLINKs .hom file obtained using Human Genome Diversity Panel. Cluster Africa.
#' PLINKs parameters were set to:
#' homozyg-window-snp 30.
#' homozyg-window-het 1.
#' homozyg-snp 30.
#' homozyg-kb 300.
#' homozyg-density 30.
#' @format A data frame with 62143 rows and 16 variables:
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
#'   \item{POP}{}
#'   \item{CLUSTER}{}
#'   \item{NAME_CLUS}{}
#' }
#' @source \url{https://cephb.fr/en/hgdp_panel.php}
"HGDP_cl_MiddleEast"

#' ROH islands in Human Genome Diversity Panel.
#'
#' This dataset contains the ROH islands of different groups of the Human Genome Diversity Panel.
#' @format A data frame with 3104 rows and 7 variables.
#' \describe{
#'   \item{Chr}{Chromosome}
#'   \item{Start}{Starting possition of the ROH island}
#'   \item{End}{Ending possition of the ROH island}
#'   \item{Length}{Length in Mb of the ROH island}
#'   \item{N_Individuals}{Number of individuals that have a ROH in that genomic location}
#'   \item{P_Individuals}{Proportion, in \%, of the individuals in the sample that have a ROH in that genomic location}
#'   \item{Population}{Name of the population}
#' }
"ROHi_HGDP"

#' Runs of Heterozygosity RHZ in Human Genome Diversity Panel.
#'
#' This dataset contains the Runs of Heterozygosity of different groups of the Human Genome Diversity Panel.
#' @format A data frame with 3104 rows and 7 variables.
#' \describe{
#'   \item{Chr}{Chromosome}
#'   \item{Start}{Starting possition of the RHZ}
#'   \item{End}{Ending possition of the RHZ}
#'   \item{Length}{Length in Mb of the RHZ}
#'   \item{N_Individuals}{Number of individuals that do not have a ROH in that genomic location, so its in RHZ}
#'   \item{P_Individuals}{Proportion, in \%, of the individuals in the sample that do not have a ROH in that genomic location}
#'   \item{Population}{Name of the population}
#' }
"RHZ_HGDP"

#' Inbreeding Depression analysis with COVID-19 data.
#'
#' This dataset contains the information needed to perform an inbreeding depression analysis with hospitalization during COVID-19.
#' @format A data frame with 7387 rows and 14 variables.
#' \describe{
#'   \item{a1hosp}{Hospitalization 0 or 1}
#'   \item{Froh}{Inbreeding genomic coefficient}
#'   \item{Fhat1}{Fgrm}
#'   \item{age}{age of the individual}
#'   \item{pc1}{Pricipal component 1}
#'   \item{pc2}{Pricipal component 2}
#'   \item{pc3}{Pricipal component 3}
#'   \item{pc4}{Pricipal component 4}
#'   \item{pc5}{Pricipal component 5}
#'   \item{pc6}{Pricipal component 6}
#'   \item{pc7}{Pricipal component 7}
#'   \item{pc8}{Pricipal component 8}
#'   \item{pc9}{Pricipal component 9}
#'   \item{pc10}{Pricipal component 10}
#' }
"COVID_ID"

#' Protein coding genes in COVID ROHi
#'
#' This dataset contains the protein coding genes found in the ROH islands of the COVID 19 dataset
#' @format A data frame with 153 rows and 7 variables.
#' \describe{
#'   \item{ensembl_gene_id}{Gene ID Ensembl format}
#'   \item{chromosome_name}{Chromosome number}
#'   \item{start_position}{start of the gene}
#'   \item{end_position}{end of the gene}
#'   \item{gene_biotype}{gene type_ protein coding gene}
#'   \item{external_gene_name}{gene name}
#'   \item{description}{description of the gene}
#' }
"COVID_genes"







