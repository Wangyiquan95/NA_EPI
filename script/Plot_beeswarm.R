# Title     : plot beeswarm
# Objective :
# Created by: yiquan
# Created on: 7/29/21
library(ggbeeswarm)
library("beeswarm")
library(ggplot2)
library(readr)
library(dplyr)
require(cowplot)
library(RColorBrewer)

textsize <- 7
level <- c('≤ 5 years','> 5 years')
colorscale  <- c(brewer.pal(8,"Accent"))
lessfive <- read_csv('result/Coevols_<5.csv')
lessfive$group <- '≤ 5 years'

greaterfive <- read_csv('result/Coevols_>5.csv')
greaterfive$group <- '> 5 years'

total <- rbind(lessfive, greaterfive)%>%
  mutate(group=factor(group, levels=level))

p <-ggplot(total, aes(group, Coevol_S,color = group)) +
  geom_beeswarm(size=2,pch=1)+theme_cowplot(12) +
  theme(plot.title=element_text(size=textsize,face="bold",hjust = 0.5),
        axis.title=element_text(size=textsize,face="bold"),
        axis.text=element_text(size=textsize,face="bold",color = 'black'),
        legend.title=element_blank(),
        legend.text=element_text(size=textsize,face="bold"),
        legend.position='none') + ylim(-1.2, 1.2)+
  scale_color_manual(values=c(colorscale[1],colorscale[2]))+
  xlab(bquote(bold(''))) +
  ylab(bquote(bold('Coevolution score')))
ggsave('graph/Coevol_byfilter.png',bg='white',p,height=4,width=4)
