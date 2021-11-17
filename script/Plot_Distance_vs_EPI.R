# Title     : Plot distance vs EPI
# Objective :
# Created by: yiquan
# Created on: 6/24/21

library(ggplot2)
library(scales)
library(RColorBrewer)
library(readr)
library(tidyr)
library(reshape)
library(stringr)
require(cowplot)
library(dplyr)
library(ggpubr)


read_epi <- function(input){
  #change the l1 to 328 and add mut1 and mut2 columnn
  epi_table  <- read_csv(input)
  epi_table <- epi_table %>%
    mutate(L1=recode(L1, `1` = '328',`2` = '329', `3` = '344',`4` = '367',`5` = '368',`6` = '369',`7` = '370'),
          L2=recode(L2, `1` = '328',`2` = '329', `3` = '344',`4` = '367',`5` = '368',`6` = '369',`7` = '370'))
  epi_table <- epi_table %>%
    mutate(mut1=paste(epi_table$AA1, epi_table$L1, sep = ''),
          mut2=paste(epi_table$AA2, epi_table$L2, sep = '')) %>%
    mutate(pair=paste(L1, L2, sep = '/'))

  mut_levels <- unique(c(epi_table$mut1,epi_table$mut2))
  epi_table  <- epi_table %>%
		    mutate(mut1=factor(mut1,levels=mut_levels)) %>%
		    mutate(mut2=factor(mut2,levels=mut_levels))
  return(epi_table)
}
PlotCor <- function (common_df,title){
	p1 <- ggplot(common_df,aes(x=distance,y=EPI)) +
      geom_point(size=1) + xlim(1,15)+ ylim(-1.6,1.6)+theme_cowplot(12) +
      theme(plot.title=element_text(size=textsize,face="bold",hjust = 0.5),
            axis.title=element_text(size=textsize,face="bold"),
            axis.text=element_text(size=textsize,face="bold",color = 'black'),
            legend.title=element_blank(),
            legend.text=element_text(size=textsize,face="bold"),
            legend.position='none') +
      ggtitle(title)+
      xlab(bquote(bold('Side-chain -- side-chain distance'))) +
      ylab(bquote(bold('Pairwise epistasis')))
	return(p1)
}

Comm_df <-function(df1,df2){
  chargelevel <- c('+-','++','--','+n','-n','nn')
  cm_df <- merge(df1,df2) %>%
    select(mut1,mut2,distance,EPI)
  cm_df$charge1 <- substr(cm_df$mut1,1,1)
  cm_df$charge2 <- substr(cm_df$mut2,1,1)
  cm_df <- cm_df %>%
    mutate(charge1=recode(charge1,'K'="+",'R'="+",'D'="-",'E'="-",.default ="n"),
          charge2=recode(charge2,'K'="+",'R'="+",'D'="-",'E'="-",.default ="n"))%>%
    mutate(charge=paste(charge1,charge2, sep = '')) %>%
    mutate(charge=recode(charge, '-+'='+-','n-'='-n','n+'='+n')) %>%
    mutate(charge=factor(charge, levels = chargelevel))
  return(cm_df)
}
#split plot
plot_split <- function (cm_df){
  cm_df1 <- cm_df %>% filter(charge=='+-')

  p1 <- PlotCor(cm_df1,'Opposite charge')
  ##++ or --
  cm_df2 <- cm_df %>% filter(charge=='++'|charge=='--')

  p2 <- PlotCor(cm_df2,'Same charge')
  ##n
  cm_df3 <- cm_df %>% filter(charge=='+n'|charge=='-n'|charge=='nn')

  p3 <- PlotCor(cm_df3,'Neutral')
  p <- list('p1'=p1,'p2'=p2,'p3'=p3)
  return(p)
}
#distance dataframe
dst_df <- read_tsv('result/SC_COG_distance.tsv') #USING sidechain center of geometry(COG)
dst_df$pair <- dst_df$pair %>%
  # str_replace_all('LYS`','')%>%
  # str_replace_all('ASN`','')%>%
  # str_replace_all('GLU`','')%>%
  # str_replace_all('SER`','')%>%
  # str_replace_all('PHE`','')%>%
  str_replace_all('-','/')
#epistasis dataframe
HK68_table <- read_epi('result/HK68_epi.csv')
Bk79_table <- read_epi('result/Bk79_epi.csv')
Bei89_table <- read_epi('result/Bei89_epi.csv')
Mos99_table <- read_epi('result/Mos99_epi.csv')
Vic11_table <- read_epi('result/Vic11_epi.csv')
HK19_table <- read_epi('result/HK19_epi.csv')
#common dataframe
common_df1 <- Comm_df(dst_df,HK68_table)
common_df2 <- Comm_df(dst_df,Bk79_table)
common_df3 <- Comm_df(dst_df,Bei89_table)
common_df4 <- Comm_df(dst_df,Mos99_table)
common_df5 <- Comm_df(dst_df,Vic11_table)
common_df6 <- Comm_df(dst_df,HK19_table)
textsize <- 7
#split them by charge
t1<-plot_split(common_df1)
t2<-plot_split(common_df2)
t3<-plot_split(common_df3)
t4<-plot_split(common_df4)
t5<-plot_split(common_df5)
t6<-plot_split(common_df6)
#multi_plot panel
p<- ggarrange(t1$p1,t1$p2,t1$p3,t3$p1,t3$p2,t3$p3,t4$p1,t4$p2,t4$p3,t5$p1,t5$p2,t5$p3,t6$p1,t6$p2,t6$p3,nrow=6,ncol=3)
ggsave('graph/Distance_vs_epi.png',bg='white',p,width = 6, height = 10)
p2<- ggarrange(t2$p1,t2$p2,t2$p3,nrow=3,ncol=1)
ggsave('graph/Distance_vs_epi_bk79.png',bg='white',p2,width = 2, height = 6)

#correlation
  cm1 <- common_df6 %>% filter(charge=='+-')
  cor(cm1$distance,cm1$EPI)
  cm2 <- common_df6 %>% filter(charge=='++'|charge=='--')
  cor(cm2$distance,cm2$EPI)
  cm3 <- common_df6 %>% filter(charge=='+n'|charge=='-n'|charge=='nn')
  cor(cm3$distance,cm3$EPI)
