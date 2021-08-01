# Title     : Plot EPI classified by charge
# Objective :
# Created by: yiquan
# Created on: 6/23/21

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

plot_heatmap <- function(epi_table,name){
  textsize <- 7

  p  <- ggplot(epi_table, aes(x=mut1,y=mut2))+
          geom_tile(aes(fill = EPI), colour = "white",size =0.1) +
          scale_fill_gradientn(name = '',colours=c("blue", "white", "red"),
                limits=c(-1.5,1.6), breaks=c(-1.5,0,1.5),
                values=rescale(c(0,1)),
                guide="colorbar",
                na.value="black") +
          scale_size_continuous(range = c(0,2.4)) +
          theme_classic() +
          theme(plot.title=element_text(size=textsize,face="bold",hjust = 0.5),
                text = element_text(size=textsize,face="bold"),
                legend.key.size = unit(0.5, 'lines'),
                legend.position = 'bottom',
                legend.direction="horizontal",
                legend.background=element_blank(),
                legend.box = "horizontal",
                axis.text.y=element_text(size=textsize,face="bold",color='black'),
                axis.title=element_blank(),
                axis.text.x=element_text(angle=90,hjust=1,vjust=0.5,size=textsize,face="bold",color='black')) +
          scale_x_discrete(position = "bottom")+ #remove legend + guides(size = FALSE,fill=FALSE)+
	  xlab("") + ggtitle(name)+
	  ylab("") +
      labs(size='')
  }

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

textsize <- 7
colorscale  <- c(brewer.pal(8,"Accent"))
HK68_table <- read_epi('result/HK68_epi.csv')
Bei89_table <- read_epi('result/Bei89_epi.csv')
Mos99_table <- read_epi('result/Mos99_epi.csv')
Vic11_table <- read_epi('result/Vic11_epi.csv')
HK19_table <- read_epi('result/HK19_epi.csv')
Bk79_table <- read_epi('result/Bk79_epi.csv')
Bk79_charge_p <- plot_heatmap(Bk79_table,'Bk79')
ggsave('graph/Bk79_epi_charge.png',Bk79_charge_p,height=2.5,width=2)

HK68_charge_p <- plot_heatmap(HK68_table,'HK68')
Bei89_charge_p <- plot_heatmap(Bei89_table,'Bei89')
Mos99_charge_p <- plot_heatmap(Mos99_table,'Mos99')
Vic11_charge_p <- plot_heatmap(Vic11_table,'Vic11')
HK19_charge_p <- plot_heatmap(HK19_table,'HK19')

charge_epi <- ggarrange(HK68_charge_p,Bk79_charge_p,Bei89_charge_p,Mos99_charge_p,Vic11_charge_p,HK19_charge_p,nrow=3,ncol=2,common.legend = TRUE,legend="bottom")
ggsave('graph/epi_charge.png',charge_epi,height=6,width=4.5)