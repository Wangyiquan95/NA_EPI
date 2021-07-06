# Title     : Plot HA/NA Mutation by year
# Objective :
# Created by: yiquan
# Created on: 6/23/21

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

mut_year_table  <- read_tsv('result/HumanH3N2_mutation_year.tsv')

textsize <- 7


p1 <- ggplot(mut_year_table,aes(x=year)) +
	geom_line(aes(y=HA,color="lightblue")) +
	geom_line(aes(y=X2,color="coral")) +
	ylim(0,90)+ scale_color_discrete(labels=c("NA", "HA"))+
	theme_cowplot(12) +
	theme(plot.title=element_text(size=textsize,face="bold",hjust = 0.5),
		  axis.title=element_text(size=textsize,face="bold"),
		  axis.text=element_text(size=textsize,face="bold",color = 'black'),
		  legend.title=element_blank(),
		  legend.text=element_text(size=textsize,face="bold"),
		  legend.position=c(0.05,0.9)) +
	xlab(bquote(bold('Year'))) +
	ylab(bquote(bold('# of amino acid mutation')))

ggsave('graph/H3N2_mutation_year.png',bg='white',p1,width = 4, height = 2.5)