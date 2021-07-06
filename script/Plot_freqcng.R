# Title     : Plot peak diagram
# Objective :
# Created by: yiquan
# Created on: 6/22/21
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

Year <- c(1968:2020)
`-344` <- c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.013513514,0.627512128,0.347738404,0.011235955,-0.008333333,0.002345309,0.005988024,0,-0.003861004,0.001782002,-0.000452643,0.000206064,0.002325581,-0.001776199,0.000567009
,0.00120919,-0.00616808,-0.00166891,0.005331622,-0.000228117,0.000741453,-0.002947178,-0.012315451,-0.010960618,-0.162850832,-0.625646657,-0.165758017)
`+369` <- c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, 0.013513514,0.64033264,0.323681936,0.02247191,-0.008333333
,-0.003642715,0.011976048,0,-0.007722008,0.007722008,0,0,0,0,0,-0.042789223,-0.833078163,-0.055689563,-0.008314204,-0.037349804,-0.015806932,0.001006612,-0.003665058,0.00149784,-0.005215912,-0.000304429,0.000960924)
textsize <- 7
colorscale  <- c(brewer.pal(12,"Set3"))

df <- data.frame(Year, `-344`, `+369`)
p1 <- ggplot(df,aes(Year,`-344`)) +
  geom_line(color=colorscale[1]) +
  theme_cowplot() +
  theme(legend.title=element_blank(),
        legend.text=element_text(size=textsize,face="bold"),
        axis.text=element_text(size=textsize,face="bold",color = 'black'),
        axis.text.x=element_text(angle=0,hjust=0.5,vjust=0.5),
        axis.title=element_text(size=textsize,face="bold"))
p2 <- ggplot(df,aes(Year, `+369`)) +
  geom_line(color=colorscale[3]) +
  theme_cowplot() +
  theme(legend.title=element_blank(),
        legend.text=element_text(size=textsize,face="bold"),
        axis.text=element_text(size=textsize,face="bold",color = 'black'),
        axis.text.x=element_text(angle=0,hjust=0.5,vjust=0.5),
        axis.title=element_text(size=textsize,face="bold"))
p <- ggarrange(p1,p2,nrow=2,ncol=1,common.legend = TRUE,legend="right")
ggsave('graph/peak_diagram.png',p,height=2,width=3)