#R code
library(ggplot2)
library(scales)
library(ggpubr)
library(RColorBrewer)
library(readr)
library(tidyr)
library(reshape)
library(stringr)
library(dplyr)
library(gridExtra)
require(cowplot)
library(ggbeeswarm)

plotaafreq <- function(aatable,poi,xlab,ylab){
  aalevels <- c('R','K','D','E','L','N','S','T','G')
  chglevels <- c('+','-','neutral')
  postable    <- aatable %>%
    filter(grepl(poi,mut)) %>%
    #sum same charge frequence
    group_by(year,charge) %>% summarise(freq = sum(freq))%>%
    #mutate(aa=factor(aa,levels=aalevels)
    mutate(charge=factor(charge,levels=chglevels))
  colorscale  <- c(brewer.pal(12,"Set3"))
  palette     <- c(colorscale[1],colorscale[3:12])
  textsize    <- 7

  p <- ggplot(postable,aes(year,freq,group=charge,color=charge)) +
    #ggplot(postable,aes(year,freq,group=aa,color=aa)) +
         geom_line() +
         scale_color_manual(values=palette,drop=FALSE) +
         theme_cowplot() +
         theme(legend.title=element_blank(),
               #legend.key.size=unit(3,"mm"),
               legend.position = 'none',
               legend.text=element_text(size=textsize,face="bold"),
               axis.text=element_text(size=textsize,face="bold"),
               axis.text.x=element_text(angle=90,hjust=0.5,vjust=0.5),
               axis.title=element_text(size=textsize,face="bold")) +
         ylab(ylab) +
         xlab(xlab) +
         scale_y_continuous(limits=c(0,1),breaks=c(0,0.5,1),labels=c('0%','50%','100%')) +
         #scale_x_continuous(breaks=c(1968,1980,1990,2000,2010,2020),labels=c(1968,1980,1990,2000,2010,2020)) +
         guides(color=guide_legend(ncol=1))
  return (p)
  }
plot_NA_antigen_charge_freq <- function (input,resi1,resi2,resi3,resi4,resi5,resi6,resi7){
  aatable <- read_tsv(input) %>% # aa level should be HumanN2Sweep_All.tsv
    mutate(aa=mapply(function(s){return(str_sub(s,-1,-1))},mut)) %>%
    mutate(mut=mapply(function(s){return(paste(str_sub(s,-1,-1),str_sub(s,2,-2),sep=''))},mut)) %>%
    mutate(charge=recode(aa,'K'="+",'R'="+",'D'="-",'E'="-",.default ="neutral"))
  p_328 <- plotaafreq(aatable,resi1,NULL,'residue 328')
  p_329 <- plotaafreq(aatable,resi2,NULL,'residue 329')
  p_344 <- plotaafreq(aatable,resi3,NULL,'residue 344')
  p_367 <- plotaafreq(aatable,resi4,NULL,'residue 367')
  p_368 <- plotaafreq(aatable,resi5,NULL,'residue 368')
  p_369 <- plotaafreq(aatable,resi6,NULL,'residue 369')
  p_370 <- plotaafreq(aatable,resi7,NULL,'residue 370')

  p <- ggarrange(p_328,p_329,p_344,p_367,p_368,p_369,p_370,nrow=7,ncol=1,legend = NULL)

  return(p)

}
plot_seasonalN1_aafreq <- function(aatable,poi,xlab,ylab){
  aalevels <- c('R','K','D','E','L','N','S','T','G')
  chglevels <- c('+','-','neutral')
  postable    <- aatable %>%
    filter(grepl(poi,mut)) %>%
    #sum same charge frequence
    group_by(year,charge) %>% summarise(freq = sum(freq))%>%
    #mutate(aa=factor(aa,levels=aalevels)
    mutate(charge=factor(charge,levels=chglevels))
  colorscale  <- c(brewer.pal(12,"Set3"))
  palette     <- c(colorscale[1],colorscale[3:12])
  textsize    <- 7

  p <- ggplot(postable,aes(year,freq,group=charge,color=charge)) +
    #ggplot(postable,aes(year,freq,group=aa,color=aa)) +
         geom_line() +
         scale_color_manual(values=palette,drop=FALSE) +
         theme_cowplot() +
         theme(legend.title=element_blank(),
               #legend.key.size=unit(3,"mm"),
               legend.position = 'none',
               legend.text=element_text(size=textsize,face="bold"),
               axis.text=element_text(size=textsize,face="bold"),
               axis.text.x=element_text(angle=90,hjust=0.5,vjust=0.5),
               axis.title=element_text(size=textsize,face="bold")) +
         ylab(ylab) +
         xlab(xlab) +
         scale_y_continuous(limits=c(0,1),breaks=c(0,0.5,1),labels=c('0%','50%','100%')) +
         scale_x_continuous(breaks=c(1940,1956,1977,2008),labels=c(1940,1956,1977,2008),limits=c(1940,2008))+
         guides(color=guide_legend(ncol=1))
  return (p)
  }
plot_seasonalN1_antigen_charge_freq <- function (input,resi1,resi2,resi3,resi4,resi5,resi6,resi7){
  aatable <- read_tsv(input) %>% # aa level should be HumanN2Sweep_All.tsv
    mutate(aa=mapply(function(s){return(str_sub(s,-1,-1))},mut)) %>%
    mutate(mut=mapply(function(s){return(paste(str_sub(s,-1,-1),str_sub(s,2,-2),sep=''))},mut)) %>%
    mutate(charge=recode(aa,'K'="+",'R'="+",'D'="-",'E'="-",.default ="neutral"))
  p_328 <- plot_seasonalN1_aafreq(aatable,resi1,NULL,'residue 328')
  p_329 <- plot_seasonalN1_aafreq(aatable,resi2,NULL,'residue 329')
  p_344 <- plot_seasonalN1_aafreq(aatable,resi3,NULL,'residue 344')
  p_367 <- plot_seasonalN1_aafreq(aatable,resi4,NULL,'residue 367')
  p_368 <- plot_seasonalN1_aafreq(aatable,resi5,NULL,'residue 368')
  p_369 <- plot_seasonalN1_aafreq(aatable,resi6,NULL,'residue 369')
  p_370 <- plot_seasonalN1_aafreq(aatable,resi7,NULL,'residue 370')

  p <- ggarrange(p_328,p_329,p_344,p_367,p_368,p_369,p_370,nrow=7,ncol=1,legend = NULL)

  return(p)

}
season_n1_p <- plot_seasonalN1_antigen_charge_freq('result/fluB/N108_Sweep_All.tsv','328','329','341','364','365','366','367')
pdm_n1_p <-plot_NA_antigen_charge_freq('result/fluB/N1pdm_Sweep_All.tsv','328','329','341','364','365','366','367')
vic_B_p <-plot_NA_antigen_charge_freq('result/fluB/B_vic_NASweep_All.tsv','328','329','344','370','371','372','373')
yma_B_p <-plot_NA_antigen_charge_freq('result/fluB/B_yam_NASweep_All.tsv','328','329','344','370','371','372','373')

p <- ggarrange(season_n1_p,pdm_n1_p,vic_B_p,yma_B_p,nrow=1,ncol=4)

ggsave('graph/fluB/B_N1_NatMutFreq_by_charge.png',p,height=7,width=5.7,bg='white')

# #output the classified position year charge frequency table
# charge_table <- read_tsv('result/N1pdm_Sweep_All.tsv') %>%
#   mutate(aa=mapply(function(s){return(str_sub(s,-1,-1))},mut)) %>%
#   mutate(mut=mapply(function(s){return(paste(str_sub(s,-1,-1),str_sub(s,2,-2),sep=''))},mut)) %>%
#   mutate(charge=recode(aa,'K'="+",'R'="+",'D'="-",'E'="-",.default ="n")) %>%
#   group_by(pos,year,charge) %>% summarise(freq = sum(freq))%>%
#   mutate(charge=paste(pos, charge, sep= ""))
# write.csv(charge_table,'result/N1pdm_NA_classified.csv')
#
# #plot classified coevolution score
# N108_Coevol  <- read_csv('result/fluB/N108_Coevols_<5.csv')
# N108_Coevol$charge1=substr(N108_Coevol$pair, 1, 1)
# N108_Coevol$charge2=substr(N108_Coevol$pair, 6, 6)
# N108_Coevol<- N108_Coevol %>%
#     mutate(paircharge=ifelse(charge1 == 'n'|charge2=='n', 'Neutral',ifelse((charge1 == '+'&charge2=='+')|(charge1 == '-'&charge2=='-'),'Same charge','Opposite charge')))
# textsize <- 7
# class_plot <- function(table,name){
#   class_plot <-ggplot(table, aes(paircharge, Coevol_S,color = paircharge)) +
#     #geom_violin() +
#     #geom_sina(size=0.2)+
#     geom_beeswarm(size=0.2,priority='random')+
#     theme_classic() +
#     theme(plot.title=element_text(size=textsize,face="bold",hjust = 0.5),
#       text = element_text(size=textsize,face="bold"),
#       legend.key.size = unit(0.5, 'lines'),
#       legend.position='none',
#       axis.text.y=element_text(size=textsize,face="bold",color='black'),
#       axis.text.x=element_text(angle=45,hjust=1,vjust=1,size=textsize,face="bold",color='black')) +
#     xlab("") +
#     ylab("Coevolution score")+ggtitle(name)+
#     labs(color='')
#   return(class_plot)
# }
#
# N108_class <- class_plot(N108_Coevol,'H1N1 seasonal strains')
#
# ggsave('graph/fluB/H1N1_seasonal_Coevol.png',N108_class,height=2,width=2,bg="white")
