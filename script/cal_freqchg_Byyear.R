# Title     : Calculate Freqchg by year for each residues
# Objective : Calculate Freqchg by year for each residues
# Created by: yiquan
# Created on: 12/20/20


library(viridis)
library(ggplot2)
library(dplyr)


#calculate the frequence change each site, each year
cal_frechg <- function(sites,t){
  pos_vec <- NULL
  year_vec <- NULL
  freqchg_vec <- NULL
  for (i in sites){
    ti <- filter(t,t$pos ==i)
    mut_num <- length(ti$year)/53
    for (row in (mut_num+1):length(ti$pos)){
      year <- ti[row,'year']
      freqchg <- abs(ti[row,'freq']-ti[row-mut_num,'freq'])
      pos_vec <-c(pos_vec,i)
      year_vec <-c(year_vec,year)
      freqchg_vec <-c(freqchg_vec,freqchg)
    }
  }
  df_freqchg <- data.frame(pos_vec,year_vec,freqchg_vec)
  return(df_freqchg)
}

#get the max frequency change of each year
cal_maxfreqchg <- function (sites,df_freqchg){
  df <- NULL
  for (i in sites){
    dfi <- df_freqchg %>%
    filter(pos_vec==i) %>%
    group_by(year_vec) %>%
    summarize(pos = mean(pos_vec),mut_freq = max(freqchg_vec))
    df <- rbind(df,dfi)
  }
  return(df)
}

#calculate the conserve residues
cal_conserve <- function (sites,df_maxfreqchg){

  conpos_vec <- NULL
  for (i in sites){
    dfi <- df_maxfreqchg %>%
      filter(pos==i) %>%
      group_by(pos) %>%
      summarize(max_frechg = max(mut_freq))
    if(dfi[1,2]<0.01){
      conpos_vec <-c(conpos_vec,i)
    }}
  return(conpos_vector)
}

t <- read.table('result/HumanN2Sweep_All.tsv', header = 1)
sites <- 1:469
df_freqchg <- cal_frechg(sites,t)
df_maxfreqchg <- cal_maxfreqchg(sites, df_freqchg)
write.csv(df_maxfreqchg,'result/HumanN2_freqchg.csv')