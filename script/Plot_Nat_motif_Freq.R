# Title     : plot natural occurring motif frequency
# Objective :
# Created by: yiquan
# Created on: 6/12/21
library(ggplot2)
library(ggridges)
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

motif_data <- read_tsv("result/Motif_ByYear.tsv")
textsize <- 9
p <- ggplot(motif_data, aes(x=year,y=motif,height=freq)) +
       geom_density_ridges(stat="identity", scale=0.9) +
       theme_ridges(grid = FALSE, center_axis_labels = TRUE) +
       theme(plot.title=element_text(size=textsize,face="bold"),
             axis.title.x=element_text(size=textsize,face="bold"),
             axis.title.y=element_blank(),
             axis.text=element_text(size=textsize,face="bold"),
             legend.title=element_blank(),
             legend.text=element_text(size=textsize,face="bold"),
             legend.position='right') +
       xlab(bquote(bold(Year))) +
       scale_x_continuous(limit=c(1968,2020))
ggsave('graph/FreqByYear.png',p,height=7,width=6)
