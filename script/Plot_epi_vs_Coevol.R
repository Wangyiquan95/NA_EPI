# Title     : Plot epi vs Coevol_S
# Objective :
# Created by: yiquan
# Created on: 6/21/21

#R code
library(ggplot2)
library(scales)
library(RColorBrewer)
library(readr)
library(tidyr)
library(reshape)
library(stringr)
library(dplyr)
library(ggrepel)
library(gridExtra)
require(cowplot)
require(reshape2)
library(ggforce)
library(ggbeeswarm)
library(ggpubr)

read_epi <- function(input){
  #change the l1 to 328 and add mut1 and mut2 columnn
  epi_table  <- read_csv(input)
  epi_table <- epi_table %>%
    mutate(L1=recode(L1, `1` = '328',`2` = '329', `3` = '344',`4` = '367',`5` = '368',`6` = '369',`7` = '370'),
          L2=recode(L2, `1` = '328',`2` = '329', `3` = '344',`4` = '367',`5` = '368',`6` = '369',`7` = '370'))
  epi_table <- epi_table %>%
    mutate(charge1=recode(AA1,'K'="+",'R'="+",'D'="-",'E'="-",.default ="n"),
          charge2=recode(AA2,'K'="+",'R'="+",'D'="-",'E'="-",.default ="n"))
  epi_table <- epi_table %>%
    mutate(mut1=paste(epi_table$charge1, epi_table$L1, sep = ''),
          mut2=paste(epi_table$charge2, epi_table$L2, sep = '')) %>%
    group_by(mut1,mut2) %>% summarise(EPI=mean(EPI)) %>%
    mutate(pair=paste(mut1, mut2, sep = '/'))

  mut_levels <- unique(c(epi_table$mut1,epi_table$mut2))
  epi_table  <- epi_table %>%
		    mutate(mut1=factor(mut1,levels=mut_levels)) %>%
		    mutate(mut2=factor(mut2,levels=mut_levels))
  return(epi_table)
}
PlotCor <- function (common_df,title){
	p1 <- ggplot(common_df,aes(x=EPI,y=Coevol_S,group=paircharge)) +
      geom_point(aes(shape=paircharge, color=paircharge),size=1) + ylim(-1.1,1.1)+ xlim(-1.2,1.6)+
      scale_shape_manual(values=c(2, 1, 5))+
      scale_color_manual(values=c(colorscale[2],colorscale[6], colorscale[5])) +theme_cowplot(12) +
      theme(plot.title=element_text(size=textsize,face="bold",hjust = 0.5),
            axis.title=element_text(size=textsize,face="bold"),
            axis.text=element_text(size=textsize,face="bold",color = 'black'),
            legend.title=element_blank(),
            legend.text=element_text(size=textsize,face="bold"),
            legend.position='right') +
      ggtitle(title)+
      xlab(bquote(bold('Pairwise epistasis'))) +
      ylab(bquote(bold('Coevolution score')))
	return(p1)
}
Comm_df <- function(df1, df2){
  cm_df <- merge(df1,df2) %>% select(pair, EPI, Coevol_S)
  cm_df$charge1=substr(cm_df$pair,1,1)
  cm_df$charge2=substr(cm_df$pair,6,6)
  cm_df <- cm_df %>%
    mutate(paircharge=ifelse(charge1 == 'n'|charge2=='n', 'Neutral',ifelse((charge1 == '+'&charge2=='+')|(charge1 == '-'&charge2=='-'),'Same charge','Opposite charge')))
  return(cm_df)
}

textsize <- 7
colorscale  <- c(brewer.pal(8,"Accent"))
Bk79_table <- read_epi('result/Bk79_epi.csv')
natural_epi  <- read_csv('result/Coevols_<5.csv')
Coevol_table <- natural_epi %>% separate(pair, c("Charge state i","Charge state j"), sep = "/")
write.csv(Coevol_table,'result/Coevol_table.csv')

common_df2 <- Comm_df(Bk79_table, natural_epi)

p2 <- PlotCor(common_df2,'Bk79')
ggsave('graph/compare_natural_epi_bk79.png',bg='white',p2,width = 4, height = 2)

#output cor of all six strain EPI with Coevol_S
HK68_table <- read_epi('result/HK68_epi.csv')
Bei89_table <- read_epi('result/Bei89_epi.csv')
Mos99_table <- read_epi('result/Mos99_epi.csv')
Vic11_table <- read_epi('result/Vic11_epi.csv')
HK19_table <- read_epi('result/HK19_epi.csv')

common_df1 <- Comm_df(HK68_table,natural_epi)
common_df3 <- Comm_df(Bei89_table,natural_epi)
common_df4 <- Comm_df(Mos99_table,natural_epi)
common_df5 <- Comm_df(Vic11_table,natural_epi)
common_df6 <- Comm_df(HK19_table,natural_epi)
p1 <- PlotCor(common_df1,'HK68')
p3 <- PlotCor(common_df3,'Bei89')
p4 <- PlotCor(common_df4,'Mos99')
p5 <- PlotCor(common_df5,'Vic11')
p6 <- PlotCor(common_df6,'HK19')
p <- ggarrange(p1,p2,p3,p4,p5,p6,nrow=3,ncol=2,common.legend = TRUE,legend="right")
ggsave('graph/compare_natural_epi.png',bg='white',p,height=7,width=6)
#pearson correlation
cor(common_df1$Coevol_S,common_df1$EPI)
cor(common_df2$Coevol_S,common_df2$EPI)
cor(common_df3$Coevol_S,common_df3$EPI)
cor(common_df4$Coevol_S,common_df4$EPI)
cor(common_df5$Coevol_S,common_df5$EPI)
cor(common_df6$Coevol_S,common_df6$EPI)
#spearman correlation
cor(common_df1$Coevol_S,common_df1$EPI,method = "spearman")
cor(common_df2$Coevol_S,common_df2$EPI,method = "spearman")
cor(common_df3$Coevol_S,common_df3$EPI,method = "spearman")
cor(common_df4$Coevol_S,common_df4$EPI,method = "spearman")
cor(common_df5$Coevol_S,common_df5$EPI,method = "spearman")
cor(common_df6$Coevol_S,common_df6$EPI,method = "spearman")