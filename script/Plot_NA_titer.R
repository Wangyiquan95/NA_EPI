# Title     : NA titer
# Objective : compare NA WT titer
# Created by: yiquan
# Created on: 6/20/21
library(readxl)
library(ggplot2)
library(dplyr)

StrainLevels <- c('HK68','Bk79','Bei89','Mos99','Vic11','HK19')
colorscale <- c(brewer.pal(8,"Set2"))
textsize <- 9
df <- read_excel('result/NA_WT_titer.xlsx')%>%
  mutate(Strain=factor(Strain,levels=StrainLevels))
# Default bar plot
p<- ggplot(df, aes(x=Strain, y=titer_mean)) +
  geom_bar(stat="identity", color="black",width=0.5) +
  scale_y_continuous(limits = c(0, 10), breaks = seq(0, 10, by = 2))+
  geom_errorbar(aes(ymin=titer_mean-titer_sd, ymax=titer_mean+titer_sd), width=.2)+
  labs(title="", x="", y = expression(bold("TCID"[50])))+
  theme_classic() +
  theme(plot.title=element_text(size=textsize,face="bold"),
        panel.grid.major = element_blank(),
        legend.position = "right",
        legend.text=element_text(size=textsize,face="bold"),
        axis.text=element_text(size=textsize,face="bold",color = 'black'),
        axis.title=element_text(size=textsize,face="bold",color='black'))
ggsave('graph/NA_titer.png',p,height=4.0,width=5)
