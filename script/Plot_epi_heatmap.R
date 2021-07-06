#R code
library(ggplot2)
library(scales)
library(RColorBrewer)
library(readr)
library(tidyr)
library(reshape)
library(stringr)
library(plyr)
library(dplyr)
library(ggrepel)
library(gridExtra)
require(cowplot)
require(reshape2)
library(ggforce)
library(ggbeeswarm)

plot_heatmap <- function(epi_table,name){
  textsize <- 7

  p  <- ggplot(epi_table, aes(x=mut1,y=mut2))+
          geom_point(aes(fill=EPI,size=abs(EPI)),color='black',pch=21) +
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
                legend.position = 'none',
                legend.direction="horizontal",
                legend.background=element_blank(),
                legend.box = "horizontal",
                axis.text.y=element_text(size=textsize,face="bold",color='black'),
                axis.title=element_blank(),
                axis.text.x.top=element_text(angle=90,hjust=1,vjust=0.5,size=textsize,face="bold",color='black')) +
          scale_x_discrete(position = "top")+ guides(size = FALSE)+#remove legend
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
    mutate(mut1=paste(epi_table$AA1, epi_table$L1, sep = ''),
          mut2=paste(epi_table$AA2, epi_table$L2, sep = ''))
  #add charge pair column
  epi_table <- epi_table %>%
    mutate(charge1=recode(AA1,'K'="+chg",'R'="+chg",'D'="-chg",'E'="-chg",.default ="neutral"),
          charge2=recode(AA2,'K'="+chg",'R'="+chg",'D'="-chg",'E'="-chg",.default ="neutral"))
  epi_table <- epi_table %>%
    mutate(charge=paste(epi_table$charge1, epi_table$charge2, sep = '/'))
  epi_table <- epi_table %>%
    mutate(charge=recode(epi_table$charge,"-chg/+chg"='+chg/-chg',"+chg/neutral"='neutral/+chg',"-chg/neutral"='neutral/-chg'))
  #filter the low value
  # epi_table$EPI[abs(epi_table$EPI) < 0.05 ] <- 0

  mut_levels <- unique(c(epi_table$mut1,epi_table$mut2))
  epi_table  <- epi_table %>%
		    mutate(mut1=factor(mut1,levels=mut_levels)) %>%
		    mutate(mut2=factor(mut2,levels=mut_levels))
  return(epi_table)
}
HK68_table <- read_epi('result/HK68_epi.csv')
HK68_p <- plot_heatmap(HK68_table,'HK68')

Bk79_table <- read_epi('result/Bk79_epi.csv')
Bk79_p <- plot_heatmap(Bk79_table,'Bk79')

Bei89_table <- read_epi('result/Bei89_epi.csv')
Bei89_p <- plot_heatmap(Bei89_table,'Bei89')

Mos99_table <- read_epi('result/Mos99_epi.csv')
Mos99_p <- plot_heatmap(Mos99_table,'Mos99')

Vic11_table <- read_epi('result/Vic11_epi.csv')
Vic11_p <- plot_heatmap(Vic11_table,'Vic11')

HK19_table <- read_epi('result/HK19_epi.csv')
HK19_p <- plot_heatmap(HK19_table,'HK19')

epi_plot <- ggarrange(HK68_p,Bk79_p,Bei89_p,Mos99_p,Vic11_p,HK19_p,nrow=3,ncol=2,common.legend = TRUE,legend="bottom")
ggsave('graph/EPI_heatmap.png',epi_plot,height=7.5,width=5)
#summary statistics by charge
print(ddply(Bk79_table,.(charge),summarise,mean_epi=mean(EPI),sd_epi=sd(EPI)))


#plot classified EPI
textsize <- 7
class_plot <- function(table,name){
  class_plot <-ggplot(table, aes(charge, EPI,color = charge)) +
    #geom_violin() +
    #geom_sina(size=0.2)+
    geom_beeswarm(size=0.2,priority='random')+
    theme_classic() +
    theme(plot.title=element_text(size=textsize,face="bold",hjust = 0.5),
      text = element_text(size=textsize,face="bold"),
      legend.key.size = unit(0.5, 'lines'),
      legend.position='none',
      axis.text.y=element_text(size=textsize,face="bold",color='black'),
      axis.text.x=element_text(angle=45,hjust=1,vjust=1,size=textsize,face="bold",color='black')) +
    xlab("") +
    ylab("Epistasis")+ggtitle(name)+
    labs(color='')
  return(class_plot)
}
HK68_class <- class_plot(HK68_table,'HK68')
Bk79_class <- class_plot(Bk79_table,'Bk79')
Bei89_class <- class_plot(Bei89_table,'Bei89')
Mos99_class <- class_plot(Mos99_table,'Mos99')
Vic11_class <- class_plot(Vic11_table,'Vic11')
HK19_class <- class_plot(HK19_table,'HK19')
class <- ggarrange(HK68_class,Bk79_class,Bei89_class,Mos99_class,Vic11_class,HK19_class,nrow=3,ncol=2,common.legend = TRUE,legend="bottom")
ggsave('graph/classified_epi.png',class,height=7.7,width=4)