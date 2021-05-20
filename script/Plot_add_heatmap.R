# Title     : Plot Additive fitness heatmap
# Objective : To plot the additive fitness heatmap with different strain background
# Created by: yiquan
# Created on: 5/20/21
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

  p  <- ggplot(epi_table, aes(x=mut,y=strains))+
          geom_point(aes(fill=add_fit,size=abs(add_fit)),color='black',pch=21) +
          scale_fill_gradientn(colours=c("blue", "white", "red"),
                limits=c(-0.7,0.7),
                values=rescale(c(0,1)),
                guide="colorbar",
                na.value="black") +
          scale_size_continuous(range = c(0,2.6)) +
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
  ggsave(output,p,height=1.8,width=3.6)
  }


add_df <- lapply(Sys.glob("result/*_add.csv"), read_csv) %>%
  bind_rows() %>%
  mutate(mut=paste(aa,pos,sep = '')) %>%
  mutate(strains=rep(c('Bei89','Bk79','HK19','HK68','Moss99','Vic11'),each=19))

mut_levels <- unique(add_df$mut)
strain_levels <- c('HK68','Bk79','Bei89','Moss99','Vic11','HK19')

add_df <- add_df %>%
  mutate(mut=factor(mut,levels=mut_levels)) %>%
  mutate(strains=factor(strains,levels=strain_levels))
plot_heatmap(add_df,'graph/add_heatmap.png')