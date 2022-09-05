##############################
#' ROH variables
#'
#'This script obtains main variables from a PLINK's home file: Sum of ROH>=1.5Mb, Number of ROH>=1.5Mb, Sum of ROH<1.5Mb,
#'Number of ROH<1.5Mb, Sum of ROH of different sizes, FROH, FoutofROH.
#' @param data A .hom file from PLINK.
#' @return A data frame with different variables.
#' @export
#'
roh_sum_id<-function(data){
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
#'This script summarize all the outcome of roh_sum_id by population.
#' @param data_1 The outcome of the function roh_sum_id.
#' @param data_2 A file with two columns: IID, pop. pop must contain each individula's population
#' @return A data frame with different variables summarize for each population.
#' @export
#'
roh_sum_pop<-function(data_1,data_2){
  mer<-merge(data_1,data_2,by="IID") |>
    dplyr::group_by(pop)|>
    dplyr::summarise(mean_Sum_long=mean(Sum_long),
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
#' Figure of the total sum of different ROH lenghts
#'
#'This script creates a figure of Sum of ROH for different length ROH classes.
#'The classification is made by continents or regions
#' @param data_1 The outcome of the function roh_sum_id.
#' @param data_2 A file with two columns: IID, pop. pop must contain each individula's population
#' @param data_3 A file with two columns: pop, cont. pop must contain all the populations present in the data_1 file, cont must contain the region of each population.
#' @return A figure of the of the total sum of ROH for different ROH size classes and continents or regions.
#' @export
#'
ROH_class_fig<-function(data_1,data_2,data_3){
  mer=merge(data_1,data_2,by="IID") |> dplyr::group_by(pop) |>
    dplyr::summarise(cl03_05=mean(cl_03_05),
                     cl05_1=mean(cl_05_1),
                     cl1_2=mean(cl_1_2),
                     cl2_4=mean(cl_2_4),
                     cl4_8=mean(cl_4_8),
                     cl8=mean(cl_8))
  df<-data.frame(mer$cl03_05,mer$cl05_1,mer$cl1_2,mer$cl2_4,mer$cl4_8,mer$cl8)
  df<-data.frame(Sum=unlist(df,use.names=FALSE))
  df<-dplyr::mutate(df,class=c(rep("cl03_05",length(unique(data_2$pop))),
                               rep("cl05_1",length(unique(data_2$pop))),
                               rep("cl1_2",length(unique(data_2$pop))),
                               rep("cl2_4",length(unique(data_2$pop))),
                               rep("cl4_8",length(unique(data_2$pop))),
                               rep("cl8",length(unique(data_2$pop)))))
  df<-dplyr::mutate(df,pop=rep(unique(KGenomes_pops$pop),6))
  df<-merge(df,data_3,by="pop")
  df<-df[order(df$cont),]
  ggplot2::ggplot(data=df, ggplot2::aes(x=class, y=Sum, group=pop,)) +
    ggplot2::geom_line(ggplot2::aes(color=cont))+
    ggplot2::geom_point(ggplot2::aes(color=cont))+
    ggplot2::theme_light()
}
##############################
#' RAW Data: Fig total Sum of ROH size classes.
#'
#'This script creates a figure of Sum of ROH for different length ROH classes.
#'The classification is made by continents or regions
#' @param data_1 The outcome of the function roh_sum_id.
#' @param data_2 A file with two columns: IID, pop. pop must contain each individual's population
#' @param data_3 A file with two columns: pop, cont. pop must contain all the populations present in the data_1 file, cont must contain the region of each population.
#' @return A data frame with the raw data needed to build the figure of the total sum of ROH for different size classes.
#' @export
#'
ROH_class_data<-function(data_1,data_2,data_3){
  mer=merge(data_1,data_2,by="IID") |> dplyr::group_by(pop) |>
    dplyr::summarise(cl03_05=mean(cl_03_05),
                     cl05_1=mean(cl_05_1),
                     cl1_2=mean(cl_1_2),
                     cl2_4=mean(cl_2_4),
                     cl4_8=mean(cl_4_8),
                     cl8=mean(cl_8))
  df<-data.frame(mer$cl03_05,mer$cl05_1,mer$cl1_2,mer$cl2_4,mer$cl4_8,mer$cl8)
  df<-data.frame(Sum=unlist(df,use.names=FALSE))
  df<-dplyr::mutate(df,class=c(rep("cl03_05",length(unique(data_2$pop))),
                               rep("cl05_1",length(unique(data_2$pop))),
                               rep("cl1_2",length(unique(data_2$pop))),
                               rep("cl2_4",length(unique(data_2$pop))),
                               rep("cl4_8",length(unique(data_2$pop))),
                               rep("cl8",length(unique(data_2$pop)))))
  df<-dplyr::mutate(df,pop=rep(unique(KGenomes_pops$pop),6))
  df<-merge(df,data_3,by="pop")
  df<-df[order(df$cont),]
  return(df)
}
##############################
#' Number vs Sum of ROH
#'
#'This script creates a figure of the Number of ROH>=1.5Mb vs. Sum of ROH>=1.5Mb.
#'It is possible to add the simulated number and sum of ROH for different consanguineous mating.
#'The dashed diagonal line represents the regression line of N vs S of ROH for two admixed populations from the 1K genomes: ACB and ASW
#' @param data_1 The outcome of the function roh_sum_id.
#' @param data_2 A file with two columns: IID, pop. pop must contain each individula's population
#' @param simul If true simulations of number and sum of ROH for differnet consanguinity matings is added.
#' @return A figure of the Number vs the Sum of ROH for ROH larger than 1.5Mb.
#' @export
#'
n_vs_sum<-function(data_1,data_2,simul=TRUE){
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
  ggplot2::ggplot(mer,ggplot2::aes(x=Sum_long,y=N_long,color=pop,shape=pop))+
    geom_abline(intercept= 4.184, slope=0.000049, linetype="dashed") +
    {if(simul)ggplot2::geom_point(data=data_sc, color="olivedrab1",shape=20,alpha = 0.25)}+
    {if(simul)ggplot2::geom_point(data=data_fc,color="yellow1",shape=20,alpha=0.25)}+
    {if(simul)ggplot2::geom_point(data=data_av,color="orangered1",shape=20,alpha=0.25)}+
    {if(simul)ggplot2::geom_point(data=data_in,color="red1",shape=20,alpha=0.25)}+
    ggplot2::geom_point()+
    ggplot2::theme_light()
}
##############################
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
#'This script searches for ROH islands in a population
#' @param POP A .hom file from PLINK with all the individuals belonging to the same group or population.
#' @param ChroNumber Chromosome number
#' @return A table with the ROH islands
#' @export
#'
get_RHOi<-function(POP,ChroNumber,population){
  nSNP=(mean(POP$NSNP[POP$KB>=1000 & POP$KB<=1100]))*0.1
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
  nsnp<-rep(nSNP,2500)
  prop<-data.p*100
  data<-data.frame(cbind(x,data.n,prop,pos,nsnp,pval))
  data.p<-subset(data,data$pval<=0.05/(lenChro/SizeWindow))
  data.re<-data.p |> dplyr::group_by(new=cumsum(c(1,diff(x)!=1))) |>
    dplyr::summarise(pos1=min(pos),pos2=max(pos),nsnp=sum(nsnp),n.ind=mean(data.n),per.ind=mean(prop))
  Chr<-rep(1,length(data.re$pos1))
  ROHi<-data.frame(cbind(Chr,data.re))
  ROHi<-dplyr::mutate(ROHi,len=(pos2-pos1)/1000000)
  ROHi<-dplyr::mutate(ROHi,pop=rep(population,length(ROHi$Chr)))
  ROHi<-dplyr::select(ROHi,Chr,pos1,pos2,len,nsnp,n.ind,per.ind)
  colnames(ROHi)<-c("Chr","Start","End","Length","N_SNP","N_Individuals","%_Individuals")
  return(ROHi)
}
################################
#' Summarize ROHi by Population
#'
#'This script summarize all the outcome of roh_sum_id by population.
#' @param mypath Path to the folder with the outcomes of the get_RHOi() function.
#' @return A data frame with different variables summarize for each population.
#' @export
#'
rohi_sum_pop<-function(mypath){
  files_list <- list.files(path=mypath, full.names=TRUE)
  dat <- data.frame()
  for (i in 1:length(files_list)) {
    dat <- rbind(dat, read.csv((files_list[i]),header=TRUE))
  }
  dat<-dat|>
    dplyr::group_by(Population)|>
    dplyr::summarise(mean_length=mean(Length),
                     sd_length=sd(Length),
                     median_Length=median(Length),
                     iqr_Length=IQR(Length),
                     max_Length=max(Length))
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
#' RUNS OF HETEROZYGOSITY i
#'
#'This script prepares the dataset to search for runs of heterozygosity in the population.
#' @param POP A .hom file from PLINK with all the individuals belonging to the same group or population.
#' @return A file ready to be used to fing RHZ
#' @export
rhc_data_org<-function(POP){
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
  return(DataFDel)
}
#library(roxygen2)
#roxygenise()




