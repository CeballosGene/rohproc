##############################
#' ROH variables
#'
#'This script extracts various valuable variables from a PLINK home file, including:
#'The total sum of ROH (Runs of Homozygosity) greater than or equal to 1.5Mb. (Sum_long).
#'The count of ROH greater than or equal to 1.5Mb (N_long).
#'The total sum of ROH less than 1.5Mb (Sum_short).
#'The count of ROH less than 1.5Mb (N_short).
#'The total sum of ROH for different size categories (0.3MB≤ROH<0.5, 0.5≤ROH<1, 1≤ROH<2, 2≤ROH4, 4≤ROH<8, ROH≥8MB).
#'The genomic inbreeding coefficient, also known as Froh.
#'The homozygosity outside of ROH, referred to as Foutroh.
#' @param data A .hom file obtained using --homozyg flag in PLINK.
#' @return A data frame with the variable previously described.
#' @export
#' @examples
#' # Obtaining basic ROH variables:
#' roh_summ_id(HGDP_hom)
#'
roh_summ_id<-function(data){
  results = data |> dplyr::group_by(IID) |>
    dplyr::summarise(Sum_long=sum(KB[KB>=1500]),
                     N_long=length(KB[KB>=1500]),
                     Sum_short=sum(KB[KB<1500]),
                     N_short=length(KB[KB<1500]),
                     cl_03_05=sum(KB[KB<=500])-sum(KB[KB<300]),
                     cl_05_1=sum(KB[KB<=1000])-sum(KB[KB<500]),
                     cl_1_2=sum(KB[KB<=2000])-sum(KB[KB<1000]),
                     cl_2_4=sum(KB[KB<=4000])-sum(KB[KB<2000]),
                     cl_4_8=sum(KB[KB<=8000])-sum(KB[KB<4000]),
                     cl_8=sum(KB[KB>8000]),
                     Froh=(sum(KB[KB>=1500]))/2881033,
                     Foutroh=(sum(KB[KB<1500]))/2881033)
  results<-as.data.frame(results)
  return(results)
}
##############################
#' Summarize by Population
#'
#'This script summarize all the outcome of ```roh_summ_id()``` by population.
#' @param data_1 The outcome of the function ```roh_summ_id()```.
#' @param data_2 A file with two columns: "IID", "pop". pop must contain each individual's population
#' @return A data frame with different variables summarize for each population.
#' @export
#' @examples
#' # Obtaining basic ROH variables:
#' rohsum<-roh_summ_id(HGDP_hom)
#' roh_summ_pop(rohsum,HGDP_pops)
#'
roh_summ_pop<-function(data_1,data_2){
  mer<-merge(data_1,data_2,by="IID") |>
    dplyr::group_by(pop)|>
    dplyr::summarise(Number=length(IID),
                     mean_Sum_long=mean(Sum_long),
                     sd_Sum_long=sd(Sum_long),
                     median_Sum_long=median(Sum_long),
                     iqr_Sum_long=IQR(Sum_long),
                     mean_N_long=mean(N_long),
                     sd_N_long=sd(N_long),
                     median_N_long=median(N_long),
                     iqr_N_long=IQR(N_long),
                     mean_Sum_short=mean(Sum_short),
                     sd_Sum_short=sd(Sum_short),
                     median_Sum_short=median(Sum_short),
                     iqr_Sum_short=IQR(Sum_short),
                     mean_N_short=mean(N_short),
                     sd_N_short=sd(N_short),
                     median_N_short=median(N_short),
                     iqr_N_short=IQR(N_short),
                     mean_Froh=mean(Froh),
                     sd_Froh=sd(Froh),
                     median_Froh=median(Froh),
                     iqr_Froh=IQR(Froh),
                     mean_Foutroh=mean(Froh),
                     sd_Foutroh=sd(Foutroh),
                     median_Foutroh=median(Foutroh),
                     iqr_Foutroh=IQR(Foutroh))
  mer<-as.data.frame(mer)
  return(mer)
}
##############################
#' Summarize by Continent
#'
#'This script summarize all the outcome of ```roh_summ_id()``` by population.
#' @param data_1 The outcome of the function ```roh_summ_id()```.
#' @param data_2 A file with two columns: "IID", "cont". cont must contain each individual's continent
#' @return A data frame with different variables summarize for each continent.
#' @export
#' @examples
#' # Obtaining basic ROH variables:
#' rohsum<-roh_summ_id(HGDP_hom)
#' pops<-merge(HGDP_pops,HGDP_cont,by="pop")
#' roh_summ_cont(rohsum,pops)
#'
roh_summ_cont<-function(data_1,data_2){
  mer<-merge(data_1,data_2,by="IID") |>
    dplyr::group_by(cont)|>
    dplyr::summarise(Number=length(IID),
                     mean_Sum_long=mean(Sum_long),
                     sd_Sum_long=sd(Sum_long),
                     median_Sum_long=median(Sum_long),
                     iqr_Sum_long=IQR(Sum_long),
                     mean_N_long=mean(N_long),
                     sd_N_long=sd(N_long),
                     median_N_long=median(N_long),
                     iqr_N_long=IQR(N_long),
                     mean_Sum_short=mean(Sum_short),
                     sd_Sum_short=sd(Sum_short),
                     median_Sum_short=median(Sum_short),
                     iqr_Sum_short=IQR(Sum_short),
                     mean_N_short=mean(N_short),
                     sd_N_short=sd(N_short),
                     median_N_short=median(N_short),
                     iqr_N_short=IQR(N_short),
                     mean_Froh=mean(Froh),
                     sd_Froh=sd(Froh),
                     median_Froh=median(Froh),
                     iqr_Froh=IQR(Froh),
                     mean_Foutroh=mean(Froh),
                     sd_Foutroh=sd(Foutroh),
                     median_Foutroh=median(Foutroh),
                     iqr_Foutroh=IQR(Foutroh))
  mer<-as.data.frame(mer)
  return(mer)
}
##############################
#' Summarize by factor
#'
#'This script summarize all the outcome of ```roh_summ_id()``` by any factor.
#' @param data_1 The outcome of the function ```roh_summ_id()```.
#' @param data_2 A file with each individual's factor
#' @param group_vars A vector with 1 or more factors do group with.
#' @return A data frame with different variables summarize for each continent.
#' @export
#' @examples
#' # Obtaining basic ROH variables:
#' rohsum<-roh_summ_id(HGDP_hom)
#' roh_summ_factor(rohsum,HGDP_pops,pop)
#'
roh_summ_factor<-function(data_1,data_2,group_vars){
  mer<-merge(data_1,data_2,by="IID") |>
    dplyr::group_by(across({{group_vars}}))|>
    dplyr::summarise(Number=length(IID),
                     mean_Sum_long=mean(Sum_long),
                     sd_Sum_long=sd(Sum_long),
                     median_Sum_long=median(Sum_long),
                     iqr_Sum_long=IQR(Sum_long),
                     mean_N_long=mean(N_long),
                     sd_N_long=sd(N_long),
                     median_N_long=median(N_long),
                     iqr_N_long=IQR(N_long),
                     mean_Sum_short=mean(Sum_short),
                     sd_Sum_short=sd(Sum_short),
                     median_Sum_short=median(Sum_short),
                     iqr_Sum_short=IQR(Sum_short),
                     mean_N_short=mean(N_short),
                     sd_N_short=sd(N_short),
                     median_N_short=median(N_short),
                     iqr_N_short=IQR(N_short),
                     mean_Froh=mean(Froh),
                     sd_Froh=sd(Froh),
                     median_Froh=median(Froh),
                     iqr_Froh=IQR(Froh),
                     mean_Foutroh=mean(Froh),
                     sd_Foutroh=sd(Foutroh),
                     median_Foutroh=median(Foutroh),
                     iqr_Foutroh=IQR(Foutroh))
  mer<-as.data.frame(mer)
  return(mer)
}
##############################
#' Figure of the total sum of different ROH lengths
#'
#'This script creates a figure of the population's average Sum of ROH for different length ROH classes.
#'(0.3MB≤ROH<0.5, 0.5≤ROH<1, 1≤ROH<2, 2≤ROH4, 4≤ROH<8, ROH≥8MB).
#'The classification is made by continents or regions.
#' @param data_1 The outcome of the function ```roh_summ_id()```.
#' @param data_2 A data frame with three columns: IID, pop and cont. pop must contain each individual's population, cont must contain each individual's region or continent.
#' @return A figure of the of the total sum of ROH for different ROH size classes and continents or regions.
#' @export
#' @examples
#' #Obtaining basic ROH variables:
#' rohsum<-roh_summ_id(HGDP_hom)
#' pops<-merge(HGDP_pops,HGDP_cont,by="pop")
#' ROH_class_fig(rohsum,pops)
#'
ROH_class_fig<-function(data_1,data_2){
  mer=merge(data_1,data_2,by="IID") |> dplyr::group_by(pop) |>
    dplyr::summarise(cl03_05=mean(cl_03_05,na.rm = TRUE),
                     cl05_1=mean(cl_05_1,na.rm = TRUE),
                     cl1_2=mean(cl_1_2,na.rm = TRUE),
                     cl2_4=mean(cl_2_4,na.rm = TRUE),
                     cl4_8=mean(cl_4_8,na.rm = TRUE),
                     cl8=mean(cl_8,na.rm = TRUE))
  df<-data.frame(mer$cl03_05,mer$cl05_1,mer$cl1_2,mer$cl2_4,mer$cl4_8,mer$cl8)
  df<-data.frame(Sum=unlist(df,use.names=FALSE))
  df<-dplyr::mutate(df,class=c(rep("cl03_05",length(unique(data_2$pop))),
                               rep("cl05_1",length(unique(data_2$pop))),
                               rep("cl1_2",length(unique(data_2$pop))),
                               rep("cl2_4",length(unique(data_2$pop))),
                               rep("cl4_8",length(unique(data_2$pop))),
                               rep("cl8",length(unique(data_2$pop)))))
  df<-dplyr::mutate(df,pop=rep(unique(mer$pop),6))
  df<-merge(df,data_2,by="pop")
  df<-df[order(df$cont),]
  ggplot2::ggplot(data=df, ggplot2::aes(x=class, y=Sum, group=pop,)) +
    ggplot2::geom_line(ggplot2::aes(color=pop))+
    ggplot2::geom_point(ggplot2::aes(color=pop))+
    ggplot2::theme_light()
}
##############################
#' Raw Data: Fig total Sum of ROH size classes.
#'
#'This script creates a table with the raw data used in the function ROH_class_fig: average Sum of ROH for different length ROH classes.
#'(0.3MB≤ROH<0.5, 0.5≤ROH<1, 1≤ROH<2, 2≤ROH4, 4≤ROH<8, ROH≥8MB).
#'The classification is made by continents or regions.
#' @param data_1 The outcome of the function ```roh_summ_id()```.
#' @param data_2 A data frame with three columns: IID, pop and cont. pop must contain each individual's population, cont must contain each individual's region or continent.
#' @return A data frame with the raw data needed to build the figure of the total sum of ROH for different size classes.
#' @export
#' @examples
#' # Obtaining basic ROH variables:
#' rohsum<-roh_summ_id(HGDP_hom)
#' pops<-merge(HGDP_pops,HGDP_cont,by="pop")
#' ROH_class_data(rohsum,pops)
#'
ROH_class_data<-function(data_1,data_2){
  mer=merge(data_1,data_2,by="IID") |> dplyr::group_by(pop) |>
    dplyr::summarise(cl03_05=mean(cl_03_05,na.rm = TRUE),
                     cl05_1=mean(cl_05_1,na.rm = TRUE),
                     cl1_2=mean(cl_1_2,na.rm = TRUE),
                     cl2_4=mean(cl_2_4,na.rm = TRUE),
                     cl4_8=mean(cl_4_8,na.rm = TRUE),
                     cl8=mean(cl_8,na.rm = TRUE))
  df<-data.frame(mer$cl03_05,mer$cl05_1,mer$cl1_2,mer$cl2_4,mer$cl4_8,mer$cl8)
  df<-data.frame(Sum=unlist(df,use.names=FALSE))
  df<-dplyr::mutate(df,class=c(rep("cl03_05",length(unique(data_2$pop))),
                               rep("cl05_1",length(unique(data_2$pop))),
                               rep("cl1_2",length(unique(data_2$pop))),
                               rep("cl2_4",length(unique(data_2$pop))),
                               rep("cl4_8",length(unique(data_2$pop))),
                               rep("cl8",length(unique(data_2$pop)))))
  df<-dplyr::mutate(df,pop=rep(unique(mer$pop),6))
  df<-merge(df,data_2,by="pop")
  df<-df[order(df$cont),]
  return(df)
}
##############################
#' Number vs Sum of ROH
#'
#'This script creates a figure of the Number of ROH>=1.5Mb vs. Sum of ROH>=1.5Mb.
#'It is possible to add the simulated number and sum of ROH for different consanguineous mating.
#'Simulated data for the number and total length of ROHs (ROH > 1.5 Mb) in the offspring resulting from various consanguineous matings are also displayed.
#'Points of varying colors represent offspring from different consanguineous relationships: green for second cousins, yellow for first cousins, orange for avuncular (including uncle-niece, aunt-nephew, and double first cousins), and red for incestuous relationships (such as brother-sister and parent-offspring).
#'Each type of consanguineous mating is depicted through 5,000 simulations.
#'The dashed diagonal line represents the regression line of N vs S of ROH for two admixed populations from the 1K genomes: ACB and ASW
#' @param data_1 The outcome of the function ```roh_sum_id()```.
#' @param data_2 A data frame with three columns: IID, pop and cont. pop must contain each individual's population, cont must contain each individual's region or continent.
#' @param color Factor variable to be introduce as the color guide: ```aes(color=,)```: pop or cont.
#' @param shape Factor variable to be introduce as the shape guide: ```aes(shape=,)```: pop or cont.
#' @param simul If true simulations of number and sum of ROH for different consanguinity mating is added.
#' @return A figure of the Number vs the Sum of ROH for ROH larger than 1.5Mb.
#' @export
#' @examples
#' # Obtaining basic ROH variables:
#' rohsum<-roh_summ_id(HGDP_hom)
#' pops<-merge(HGDP_pops,HGDP_cont,by="pop")
#' n_vs_sum(rohsum,pops,pop,cont,simul=TRUE)
#'
n_vs_sum<-function(data_1,data_2,color,shape,simul=TRUE){
  Sum_long<-rnorm(5000,0.0152,0.009)*2881033
  N_long<-rnorm(5000,4.81,2.5)
  data_sc<-as.data.frame(cbind(Sum_long,N_long));data_sc[data_sc < 0] <- 0
  Sum_long<-rnorm(5000,0.0625,0.011)*2881033
  N_long<-rnorm(5000,14.8,2.5)
  data_fc<-as.data.frame(cbind(Sum_long,N_long))
  Sum_long<-rnorm(5000,0.1252,0.02)*2881033
  N_long<-rnorm(5000,20.5,4.5)
  data_av<-as.data.frame(cbind(Sum_long,N_long))
  Sum_long<-rnorm(5000,0.25,0.03)*2881033
  N_long<-rnorm(5000,40.7,4.5)
  data_in<-as.data.frame(cbind(Sum_long,N_long))
  mer<-merge(data_1,data_2,by="IID")
  ggplot2::ggplot(mer,ggplot2::aes(x=Sum_long,y=N_long,color={{color}},shape={{shape}}))+
    ggplot2::geom_abline(intercept= 0.4866, slope=0.000438, linetype="dashed") +
    {if(simul)ggplot2::geom_point(data=data_sc, color="olivedrab1",shape=20,alpha = 0.25)}+
    {if(simul)ggplot2::geom_point(data=data_fc,color="yellow1",shape=20,alpha=0.25)}+
    {if(simul)ggplot2::geom_point(data=data_av,color="orangered1",shape=20,alpha=0.25)}+
    {if(simul)ggplot2::geom_point(data=data_in,color="red1",shape=20,alpha=0.25)}+
    ggplot2::geom_point()+
    ggplot2::theme_light()
}
##############################
# FIS vs FROH
#'
#'This script creates a scatter plot of FIS vs FROH.
#' @param data_1 The outcome of the function ```roh_sum_id()```.
#' @param data_2 A data frame with three columns: IID, pop and cont. pop must contain each individual's population, cont must contain each individual's region or continent.
#' @param color Factor variable to be introduce as the color guide: ```aes(color=,)```: pop or cont.
#' @param shape Factor variable to be introduce as the shape guide: ```aes(shape=,)```: pop or cont.
#' @param simul If true simulations of number and sum of ROH for different consanguinity mating is added.
#' @return A figure of the Number vs the Sum of ROH for ROH larger than 1.5Mb.
#' @export
#' @examples
#' # Obtaining basic ROH variables:
#' rohsum<-roh_summ_id(HGDP_hom)
#' rohsum<-merge(rohsum,HGDP_het,by="IID")
#' pops<-merge(HGDP_pops,HGDP_cont,by="pop")
#' fis_vs_froh(rohsum,pops,pop,cont)
#'
fis_vs_froh<-function(data_1,data_2,color,shape){
  mer<-merge(data_1,data_2,by="IID")
  ggplot2::ggplot(mer, ggplot2::aes(x=Froh, y=Fis, color={{color}},shape={{shape}})) +
    ggplot2::geom_point()+
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed",color="red") +
    ggplot2::geom_abline(slope=1,intercept = 0, linetype = "dashed",color="red")+
    ggplot2::theme_light()
}
############################
# FUNCTIONS ROHi ----
roh_island<-function(pop,chr,p1,p2){
  names(pop)<-tolower(names(pop))
  a<-pop[pop$chr==chr,]
  island<-subset(a,pos1<=p1 & pos2>=p2)
  n<-length(unique(island$iid))/length(unique(pop$iid))
  return(n)
}
poisson.roh_island<-function(pop,chr,p1,p2){
  names(pop)<-tolower(names(pop))
  a<-pop[pop$chr==chr,]
  island<-subset(a,pos1<=p1 & pos2>=p2)
  n<-length(unique(island$iid))
  return(n)
}
##############################
#' Islands of ROH
#'
#'This script searches for ROH islands in a population. The script may take some time to finish the 22 Chr.
#' @param POP A .hom file from PLINK with all the individuals belonging to the same group or population.
#' @param ChroNumber Chromosome number.
#' @param population Name of the population in quotation marks ("").
#' @return A table with the ROH islands.
#' @export
#'
get_RHOi<-function(POP,ChroNumber,population){
  SizeWindow=10000
  if(ChroNumber==1){lenChro=250000000}
  if(ChroNumber==2){lenChro=250000000}
  if(ChroNumber==3){lenChro=200000000}
  if(ChroNumber==4){lenChro=191000000}
  if(ChroNumber==5){lenChro=182000000}
  if(ChroNumber==6){lenChro=171000000}
  if(ChroNumber==7){lenChro=160000000}
  if(ChroNumber==8){lenChro=146000000}
  if(ChroNumber==9){lenChro=139000000}
  if(ChroNumber==10){lenChro=133900000}
  if(ChroNumber==11){lenChro=136000000}
  if(ChroNumber==12){lenChro=134000000}
  if(ChroNumber==13){lenChro=115000000}
  if(ChroNumber==14){lenChro=108000000}
  if(ChroNumber==15){lenChro=102000000}
  if(ChroNumber==16){lenChro=91000000}
  if(ChroNumber==17){lenChro=84000000}
  if(ChroNumber==18){lenChro=81000000}
  if(ChroNumber==19){lenChro=59000000}
  if(ChroNumber==20){lenChro=64000000}
  if(ChroNumber==21){lenChro=49000000}
  if(ChroNumber==22){lenChro=52000000}
  data.n = apply(data.frame(seq(0,lenChro-SizeWindow,SizeWindow), seq(SizeWindow,lenChro,SizeWindow)), MARGIN=1,
                 function(x,y,z,a) poisson.roh_island(POP, ChroNumber, x[1], x[2]))
  data.p = apply(data.frame(seq(0,lenChro-SizeWindow,SizeWindow), seq(SizeWindow,lenChro,SizeWindow)), MARGIN=1,
                 function(x,y,z,a) roh_island(POP, ChroNumber, x[1], x[2]))
  av<-mean(data.n)
  pval<-ppois(data.n,lambda=av,lower=FALSE)
  x<-seq(1:round(lenChro/SizeWindow))
  pos<-(x*SizeWindow)
  prop<-data.p*100
  data<-data.frame(cbind(x,data.n,prop,pos,pval))
  data.p<-subset(data,data$pval<=0.05/(lenChro/SizeWindow))
  data.re<-data.p |> dplyr::group_by(new=cumsum(c(1,diff(x)!=1))) |>
    dplyr::summarise(pos1=min(pos),pos2=max(pos),n.ind=mean(data.n),per.ind=mean(prop))
  Chr<-rep(ChroNumber,length(data.re$pos1))
  ROHi<-data.frame(cbind(Chr,data.re))
  ROHi<-dplyr::mutate(ROHi,len=(pos2-pos1)/1000000)
  ROHi<-dplyr::mutate(ROHi,pop=rep(population,length(ROHi$Chr)))
  ROHi<-dplyr::select(ROHi,Chr,pos1,pos2,len,n.ind,per.ind,pop)
  colnames(ROHi)<-c("Chr","Start","End","Length","N_Individuals","P_Individuals","Population")
  return(ROHi)
}
################################
#' Summarize ROHi by Population
#'
#'This script summarize all the outcome of get_ROHi by population.
#' @param mypath Path to the folder with the outcomes of the get_RHOi() function.
#' @return A data frame with different variables summarize for each population.
#' @export
#'
rohi_summ_pop<-function(mypath){
  files_list <- list.files(path=mypath, full.names=TRUE)
  dat <- data.frame()
  for (i in 1:length(files_list)) {
    dat <- rbind(dat, read.csv((files_list[i]),header=TRUE))
  }
  dat<-dat|>
    dplyr::group_by(Population)|>
    dplyr::summarise(Number=length(Chr),
                     mean_length=mean(Length),
                     sd_length=sd(Length),
                     median_Length=median(Length),
                     iqr_Length=IQR(Length),
                     max_Length=max(Length),
                     mean_N_Individuals=mean(N_Individuals),
                     sd_N_Individuals=sd(N_Individuals),
                     median_N_Individuals=median(N_Individuals),
                     iqr_N_Individuals=IQR(N_Individuals),
                     max_N_Individuals=max(N_Individuals),
                     mean_P_Individuals=mean(P_Individuals),
                     sd_P_Individuals=sd(P_Individuals),
                     median_P_Individuals=median(P_Individuals),
                     iqr_P_Individuals=IQR(P_Individuals),
                     max_P_Individuals=max(P_Individuals))
  out<-as.data.frame(dat)
  return(out)
}
##############################
# FUNCTIONS RHZ ----
RemoveBlackList<-function( Start,End, Chro, IID, blacklistChro){
  Case1<-(blacklistChro$Start>Start & blacklistChro$Start<End) | (blacklistChro$Stop>Start & blacklistChro$Stop<End)
  Case2<-Start>blacklistChro$Start & End<blacklistChro$Stop
  if(any(Case2))return(data.frame(Chro=c(), IID=c(),V1=c(), V2=c()))
  if(any(Case1)){
    blacklistChroCut<-blacklistChro[Case1,]
    blacklistChroCut<-blacklistChroCut[order(blacklistChroCut$Start),]
    Res<-data.frame(Chro=Chro, IID=IID,V1=c(Start,blacklistChroCut$Stop), V2=c(blacklistChroCut$Start,End))
    if(Start>blacklistChroCut$Start[1])Res<-Res[-1,]
    if(End<blacklistChroCut$Stop[length(blacklistChroCut$Stop)])Res<-Res[-nrow(Res),]
    return(Res)
  }else{
    return(data.frame(Chro=Chro, IID=IID,V1=Start, V2=End))
  }
}
##############################
#' RUNS OF HETEROZYGOSITY I
#'
#'This script prepares the dataset to search for runs of heterozygosity in the population.
#' @param POP A .hom file from PLINK with all the individuals belonging to the same group or population.
#' @return A file ready to be used to fing RHZ.
#' @export
rhc_data_org<-function(POP){
  positions=help_RHZ[[1]]
  blacklist=help_RHZ[[2]]
  POP<- POP[with(POP,order(CHR,IID,POS1)),]
  Cmt=1
  for(Chro in 1:22){
    for(Ind in unique(POP$IID)){
      DataChro<-data.frame(X1=c(),X2=c())
      Balise<-POP$CHR==Chro & POP$IID==Ind
      if(any(Balise)){
        if(positions[positions$chr==Chro,"snp.P1"]!=POP[Balise,"POS1"][1]){
          Begin<-cbind(positions[positions$chr==Chro,"snp.P1"], POP[Balise,"POS1"][1]-1)
          DataChro<-rbind(DataChro,Begin)
        }
        if(nrow(POP[Balise,])>1){
          DataChro2<-cbind(POP[Balise,"POS2"][1:c(length(POP[Balise,"POS2"])-1)]+1,POP[Balise,"POS1"][2:length(POP[Balise,"POS1"])]-1)
          DataChro<-rbind(DataChro,DataChro2)
        }
        if(positions[positions$chr==Chro,"snp.P2"]!=POP[Balise,"POS2"][length(POP[Balise,"POS2"])]){
          End<-cbind( POP[Balise,"POS2"][length(POP[Balise,"POS2"])]+1, positions[positions$chr==Chro,"snp.P2"])
          DataChro<-rbind(DataChro,End)
        }
        PosCentro<-which(DataChro[,1]<positions[positions$chr==Chro,"cen.pos1"] & DataChro[,2]>positions[positions$chr==Chro,"cen.pos2"])
        if(length(PosCentro)==1){
          DataCentro<-rbind(cbind(DataChro[PosCentro,1], positions[positions$chr==Chro,"cen.pos1"]),cbind(positions[positions$chr==Chro,"cen.pos2"],DataChro[PosCentro,2]))
          DataChro<-DataChro[-PosCentro,]
          DataChro<-rbind(DataChro,DataCentro)
          DataChro<-DataChro[order(DataChro[,1]),]
        }
        if(length(PosCentro)>1){
          cat("Error :", Chro)
        }
        DataChro<-data.frame(Chro=Chro,IID=Ind,DataChro)
        if(Cmt==1)DataF<-DataChro
        else DataF<-rbind(DataF,DataChro)
        Cmt=Cmt+1
      }
    }
  }
  Cmt<-1
  for(Chro in 1:22){
    blacklistChro<-blacklist[blacklist$Chr==Chro,]
    DataFChro<-DataF[DataF$Chro==Chro,]
    for(eachRoh in 1:nrow(DataFChro)){
      DataTemp<-RemoveBlackList(DataFChro$V1[eachRoh], DataFChro$V2[eachRoh], DataFChro$Chro[eachRoh],DataFChro$IID[eachRoh],blacklistChro)
      if(Cmt==1)DataFDel<-DataTemp
      else DataFDel<-rbind(DataFDel,DataTemp)
      Cmt<-Cmt+1
    }
  }
  colnames(DataFDel) <- c("CHR","IID","POS1","POS2")
  cluster<-DataFDel
  cluster$KB<-cluster$POS2 - cluster$POS1
  return(cluster)
}
##############################
#' RUNS OF HETEROZYGOSITY II
#'
#'This script searches for Regions of Heterozygosity in a population
#' @param DF The data frame obtained from the function ```rhc_data_org()```.
#' @param ChroNumber Chromosome number
#' @param PR Percentage of population without ROH in a particular window. 100%: no individual with homozygous haplotype
#' in that window. 90%: 10% of people with homozygous regions in that window.
#' @param population Name of the population in quotation marks ("").
#' @return A table with the RHZ.
#' @export
#'
get_RHZ<-function(DF,ChroNumber,PR,population){
  SizeWindow=10000
  if(ChroNumber==1){lenChro=250000000}
  if(ChroNumber==2){lenChro=250000000}
  if(ChroNumber==3){lenChro=200000000}
  if(ChroNumber==4){lenChro=191000000}
  if(ChroNumber==5){lenChro=182000000}
  if(ChroNumber==6){lenChro=171000000}
  if(ChroNumber==7){lenChro=160000000}
  if(ChroNumber==8){lenChro=146000000}
  if(ChroNumber==9){lenChro=139000000}
  if(ChroNumber==10){lenChro=133900000}
  if(ChroNumber==11){lenChro=136000000}
  if(ChroNumber==12){lenChro=134000000}
  if(ChroNumber==13){lenChro=115000000}
  if(ChroNumber==14){lenChro=108000000}
  if(ChroNumber==15){lenChro=102000000}
  if(ChroNumber==16){lenChro=91000000}
  if(ChroNumber==17){lenChro=84000000}
  if(ChroNumber==18){lenChro=81000000}
  if(ChroNumber==19){lenChro=59000000}
  if(ChroNumber==20){lenChro=64000000}
  if(ChroNumber==21){lenChro=49000000}
  if(ChroNumber==22){lenChro=52000000}
  data.n = apply(data.frame(seq(0,lenChro-SizeWindow,SizeWindow), seq(SizeWindow,lenChro,SizeWindow)), MARGIN=1,
                 function(x,y,z,a) poisson.roh_island(DF, ChroNumber, x[1], x[2]))
  data.p = apply(data.frame(seq(0,lenChro-SizeWindow,SizeWindow), seq(SizeWindow,lenChro,SizeWindow)), MARGIN=1,
                 function(x,y,z,a) roh_island(DF, ChroNumber, x[1], x[2]))
  x<-seq(1:round(lenChro/SizeWindow))
  pos<-(x*SizeWindow)
  prop<-data.p*100
  data<-data.frame(cbind(x,data.n,prop,pos))
  data.p<-subset(data,data$prop>=PR)
  data.re<-data.p |> dplyr::group_by(new=cumsum(c(1,diff(x)!=1))) |>
    dplyr::summarise(pos1=min(pos),pos2=max(pos),n.ind=mean(data.n),per.ind=mean(prop))
  Chr<-rep(ChroNumber,length(data.re$pos1))
  OUT<-data.frame(cbind(Chr,data.re))
  OUT<-dplyr::mutate(OUT,len=(pos2-pos1)/1000000)
  OUT<-dplyr::mutate(OUT,pop=rep(population,length(OUT$Chr)))
  OUT<-dplyr::select(OUT,Chr,pos1,pos2,len,n.ind,per.ind,pop)
  colnames(OUT)<-c("Chr","Start","End","Length","N_Individuals","P_Individuals","Population")
  return(OUT)
}
################################
#' Summarize RHZ by Population
#'
#'This script summarize all the outcome of RHZ by population.
#' @param mypath Path to the folder with the outcomes of the ```get_RHZ()``` function.
#' @return A data frame with different variables summarize for each population.
#' @export
#'
rhz_summ_pop<-function(mypath){
  files_list <- list.files(path=mypath, full.names=TRUE)
  dat <- data.frame()
  for (i in 1:length(files_list)) {
    dat <- rbind(dat, read.csv((files_list[i]),header=TRUE))
  }
  dat<-dat|>
    dplyr::group_by(Population)|>
    dplyr::summarise(Number=length(Chr),
                     mean_length=mean(Length),
                     sd_length=sd(Length),
                     median_Length=median(Length),
                     iqr_Length=IQR(Length),
                     max_Length=max(Length),
                     mean_N_Individuals=mean(N_Individuals),
                     sd_N_Individuals=sd(N_Individuals),
                     median_N_Individuals=median(N_Individuals),
                     iqr_N_Individuals=IQR(N_Individuals),
                     max_N_Individuals=max(N_Individuals),
                     mean_P_Individuals=mean(P_Individuals),
                     sd_P_Individuals=sd(P_Individuals),
                     median_P_Individuals=median(P_Individuals),
                     iqr_P_Individuals=IQR(P_Individuals),
                     max_P_Individuals=max(P_Individuals))
  out<-as.data.frame(dat)
  return(out)
}
################################
#' Searching for Protein Coding Genes
#'
#'This script harvest protein coding genes from the ROHi or RHZ regions. It uses biomart package
#' @param DATA A file obtained from the functions ```get_ROHi()``` or ```get_RHZ()```.
#' @return A table with the Protein coding genes.
#' @export
get_Prot<-function(DATA){
  mart <- biomaRt::useMart("ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl", host="www.ensembl.org")
  attributes = biomaRt::listAttributes(mart)
  attributes <- c("ensembl_gene_id","chromosome_name","start_position","end_position",
                  "gene_biotype","external_gene_name","description")
  filters <- c("chromosome_name","start","end")
  ann<-list()
  gen<-list()
  for (i in seq(1,length(DATA$Chr),1)){
    ann[[i]]<-list(chromosome=DATA[i,1],start=DATA[i,2],end=DATA[i,3])
    gen[[i]]<-biomaRt::getBM(attributes=attributes, filters=filters, values=ann[[i]], mart=mart)
    gen[[i]]<-gen[[i]][gen[[i]]$gene_biotype=="protein_coding",]
  }
  genes<- plyr::ldply(gen, data.frame)
  return(genes)
}
##############################
#' Common and Unique regions
#'
#'This function finds the unique and common ROH islands or regions of heterozygosity between two
#'populations.
#' @param data_1 A file obtained from the functions ```get_ROHi()``` or ```get_RHZ()``` for the first population.
#' @param data_2 A file obtained from the functions ```get_ROHi()``` or ```get_RHZ()``` for the second population.
#' @param pop target population you want to get the outcome, in quotation marks ("").
#' @param class Type of outcomes we want. Two choices: "Unique" or "Common". n quotation marks ("").
#' @return Pairwise common and Unique regions
#' @export
comm_uni<-function(data_1,data_2,pop,class){
  Data<-rbind(data_1,data_2)
  CmtChro=1
  Res=list()
  Cmt1<-1
  for(Chro in unique(Data$Chr)){
    DataChro=Data[Data$Chr==Chro,]
    rownames(DataChro) <- seq_len(nrow(DataChro))
    ints <- tapply(1:nrow(DataChro),DataChro$Population,function(v)
      intervals::Intervals(as.matrix(DataChro[v,c("Start","End")]),
                           closed=c(FALSE,FALSE),
                           type="Z"))
    pops <- unique(DataChro$Population)
    popidx <- lapply(pops,function(x) which(DataChro$Population==x))
    names(popidx) <- pops
    sets <- expand.grid(pops,pops,stringsAsFactors = FALSE)
    sets <- sets[sets$Var1!=sets$Var2,]
    olap <- lapply(1:nrow(sets),function(i)
      intervals::interval_overlap(ints[[sets$Var1[i]]],ints[[sets$Var2[i]]]))
    olap <- lapply(1:nrow(sets),function(i) {
      df<-as.data.frame(olap[[i]],stringsAsFactors=FALSE)
      df$Start <- as.numeric(rownames(df))
      df$End <- sapply(1:nrow(df),function(j) popidx[[sets$Var2[i]]][df[j,1][[1]][1]])
      return(df)})
    olap <- do.call(rbind,olap)[,-1]
    olap$olaps <- !is.na(olap$End)
    olap <- data.frame(minoverlap=tapply(olap$olaps,olap$Start,min),maxoverlap=tapply(olap$olaps,olap$Start,max))
    olap$rowno <- as.numeric(rownames(olap))
    uniques <- DataChro[olap$rowno[olap$maxoverlap==0],]
    commons <- DataChro[olap$rowno[olap$minoverlap>0],]
    Res[[Chro]]=list('uniques'=uniques, 'commons'=commons)
    if(Cmt1==1){
      CommonAll<-commons
      UniqueAll<-uniques
    }else{
      CommonAll<-rbind(CommonAll,commons)
      UniqueAll<-rbind(UniqueAll ,uniques)
    }
    Cmt1<-Cmt1+1
  }
  CommonAll<-dplyr::mutate(CommonAll,type=rep("Common",length(CommonAll$Chr)))
  UniqueAll<-dplyr::mutate(UniqueAll,type=rep("Unique",length(UniqueAll$Chr)))
  com_uni<-rbind(CommonAll,UniqueAll)
  res<-subset(com_uni,com_uni$type==class&com_uni$Population==pop)
  return(res)
}
##############################
#' ROHi Representation
#'
#'This script creates a genomic representation of the ROH islands or Regions of Heterozygosity (RHZ).
#' @param data_1 A file obtained from the functions ```get_ROHi()``` or ```get_RHZ()``` for a population
#' @param pop target population you want to get the outcome, in quotation marks ("").
#' @return Figure of the ROHi or RHZ in a popualtion
#' @export
rohi_genome<-function(data_1,pop){
  dat_1<-subset(data_1,data_1$Population=="pop")
  positions=help_RHZ[[1]]
  d.chr=help_RHZ[[3]]
  ggplot2::ggplot(data_1)+
    ggplot2::geom_segment(ggplot2::aes(y = Chr, yend = Chr, x = 0, xend = long), data = d.chr) +
    ggplot2::geom_segment(ggplot2::aes(y = chr, yend = chr, x = cen.pos1, xend = cen.pos2, color='Centromer'),data=positions, lwd = 3)+
    ggplot2::geom_segment(ggplot2::aes(y = Chr, yend = Chr, x = Start, xend = End, color='ROHi'), lwd = 3)+
    ggplot2::scale_color_manual(name='',values=c('ROHi'='red4',"Centromer"="black"))+
    ggplot2::scale_y_discrete( limits=c(1:22)) +
    ggplot2::theme_light()
}
##############################
#' ROHi + ROHZ Representation
#'
#'This script creates a genomic representation of ROH islands and regions of heterozygosity (RHZ).
#' @param data_1 A file obtained from the functions ```get_ROHi()```.
#' @param data_2 A file obtained from the functions ```get_RHZ()```.
#' @param pop population to be represented.
#' @return Figure of the ROHi and RHZ in a popualtion
#' @export
rohi_rohz_genome<-function(data_1,data_2,pop){
  dat_1<-subset(data_1,data_1$Population=="pop")
  dat_2<-subset(data_2,data_2$Population=="pop")
  positions=help_RHZ[[1]]
  d.chr=help_RHZ[[3]]
  ggplot2::ggplot(data_1)+
    ggplot2::geom_segment(ggplot2::aes(y = Chr, yend = Chr, x = 0, xend = long), data = d.chr) +
    ggplot2::geom_segment(ggplot2::aes(y = Chr, yend = Chr, x = Start, xend = End, color='ROHi'), lwd = 3)+
    ggplot2::geom_segment(ggplot2::aes(y = chr, yend = chr, x = cen.pos1, xend = cen.pos2, color='Centromer'),data=positions, lwd = 3)+
    ggplot2::geom_segment(ggplot2::aes(y = Chr, yend = Chr, x = Start, xend =End, color = 'RHZ'),data=data_2, lwd = 3)+
    ggplot2::scale_color_manual(name='',values=c('ROHi'='red3','Centromer'='black','RHZ'='springgreen4'))+
    ggplot2::labs(y= "Chromosome", x = "Length of the Chromosome")+
    ggplot2::theme_light()
}
##############################
#library(roxygen2)
#roxygenise()

