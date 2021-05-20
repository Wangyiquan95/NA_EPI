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
epi_table  <- read_csv('result/HK19_epi.csv')
epi_table <- epi_table %>%
  mutate(L1=recode(L1, `1` = '328',`2` = '329', `3` = '344',`4` = '367',`5` = '368',`6` = '369',`7` = '370'),
         L2=recode(L2, `1` = '328',`2` = '329', `3` = '344',`4` = '367',`5` = '368',`6` = '369',`7` = '370'))
epi_table <- epi_table %>%
  mutate(mut1=paste(epi_table$AA1, epi_table$L1, sep = ''),
         mut2=paste(epi_table$AA2, epi_table$L2, sep = ''))
#filter the low value
epi_table$EPI[abs(epi_table$EPI) < 0.05 ] <- 0

mut_levels <- unique(c(epi_table$mut1,epi_table$mut2))
epi_table  <- epi_table %>%
		  mutate(mut1=factor(mut1,levels=mut_levels)) %>%
		  mutate(mut2=factor(mut2,levels=mut_levels))
plot_heatmap(epi_table,'graph/HK19_epi_heatmap.png')

