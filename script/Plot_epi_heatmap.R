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

plot_heatmap <- function(epi_table,output){
  textsize <- 7

  p  <- ggplot(epi_table, aes(x=mut1,y=mut2))+
          geom_point(aes(fill=EPI,size=abs(EPI)),color='black',pch=21) +
          scale_fill_gradientn(colours=c("blue", "white", "red"),
                limits=c(-1.5,1.6),
                values=rescale(c(0,1)),
                guide="colorbar",
                na.value="black") +
          scale_size_continuous(range = c(0,2.4)) +
          theme_classic() +
          theme(#panel.border = element_rect(colour = "black", fill=NA, size=1.2),
                text = element_text(size=textsize,face="bold"),
                legend.key.size = unit(0.5, 'lines'),
                axis.text.y=element_text(size=textsize,face="bold",color='black'),
                axis.title=element_blank(),
                axis.text.x.top=element_text(angle=90,hjust=1,vjust=0.5,size=textsize,face="bold",color='black')) +
          scale_x_discrete(position = "top") +
	  xlab("") + 
	  ylab("")
  ggsave(output,p,height=2.0,width=2.8)
  }
#change the l1 to 328 and add mut1 and mut2 columnn
epi_table  <- read_csv('result/Vic11_epi.csv')
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
plot_heatmap(epi_table,'graph/Vic11_epi_heatmap.png')

#plot classified EPI
textsize <- 7
class_plot <-ggplot(epi_table, aes(charge, EPI)) +
  #geom_violin() +
  #geom_sina(size=0.2)+
  geom_beeswarm(size=0.2,priority='random')+
  theme_classic() +
  theme(
    text = element_text(size=textsize,face="bold"),
    legend.key.size = unit(0.5, 'lines'),
    axis.text.y=element_text(size=textsize,face="bold",color='black'),
    axis.title=element_blank(),
    axis.text.x=element_text(angle=45,hjust=1,vjust=1,size=textsize,face="bold",color='black')) +
  xlab("") +
  ylab("Epistasis")
ggsave('graph/Vic11_classified_epi.png',class_plot,height=2.0,width=2.8)