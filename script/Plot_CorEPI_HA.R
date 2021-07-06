# Title     : correlate EPI among diferent strains on HA
# Objective :
# Created by: yiquan
# Created on: 6/7/21

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
library(data.table)
library(GGally)
library(e1071)
library(sinaplot)
library(ggforce)
require(cowplot)



multmerge <- function(mypath,key){
  df <- list.files(path =mypath,pattern = key,full.names = T) %>%
    lapply(read_csv) %>% bind_rows
  return(df)
}

###HA model EPI
HK68_epi <- read.csv('HA_result/HA_HK68_epi.csv') %>% mutate(HK68=EPI)
Bk79_epi <- read.csv('HA_result/HA_Bk79_epi.csv')%>% mutate(Bk79=EPI)
Bei89_epi <- read.csv('HA_result/HA_Bei89_epi.csv')%>% mutate(Bei89=EPI)
Mos99_epi <- read.csv('HA_result/HA_Mos99_epi.csv')%>% mutate(Mos99=EPI)
Bris07_epi <- read.csv('HA_result/HA_Bris07_epi.csv')%>% mutate(Bris07=EPI)
ND16_epi <- read.csv('HA_result/HA_NDako16_epi.csv')%>% mutate(ND16=EPI)
#merge df
com_epi_df <- bind_cols(list(HK68_epi,Bk79_epi,Bei89_epi,Mos99_epi,Bris07_epi,ND16_epi)) %>%
  select(2:5,7,14,21,28,35,42)
textsize <- 7


p  <- ggpairs(com_epi_df,columns=5:10,lower =list(continuous=wrap(ggally_points,size=0.1,alpha=1)),
		upper = list(continuous="blank"),
		diag  = list(continuous ="blank")) +
  theme_cowplot(9) +
  theme(plot.title=element_text(size=textsize,face="bold"),
            panel.grid.major = element_blank(),
            legend.text=element_text(size=textsize,face="bold"),
            axis.title=element_text(size=textsize,face="bold"),
            axis.text=element_text(size=textsize,face="bold"))+
  theme(plot.title=element_text(size=textsize,face="bold"),
            panel.grid.major = element_blank(),
            legend.text=element_text(size=textsize,face="bold"),
            axis.title=element_text(size=textsize,face="bold"),
            axis.text=element_text(size=textsize,face="bold")) +
  scale_x_continuous(n.breaks = 3.5) +
  scale_y_continuous(n.breaks = 4)
ggsave('graph/HA_compare_EPI.png',p,width = 4, height = 4)
cor_epi <- cor(data.table(select(com_epi_df, 5:10)))
print (cor(data.table(select(com_epi_df, 5:10))))

#####     compare additive correlation      ##########

###Load HA model add
HK68_add <- read.csv('HA_result/HA_HK68_add.csv') %>% mutate(HK68=add_fit)
Bk79_add <- read.csv('HA_result/HA_Bk79_add.csv')%>% mutate(Bk79=add_fit)
Bei89_add <- read.csv('HA_result/HA_Bei89_add.csv')%>% mutate(Bei89=add_fit)
Mos99_add <- read.csv('HA_result/HA_Mos99_add.csv')%>% mutate(Mos99=add_fit)
Bris07_add <- read.csv('HA_result/HA_bRIS07_add.csv')%>% mutate(Bris07=add_fit)
ND16_add <- read.csv('HA_result/HA_ndAKO16_add.csv')%>% mutate(ND16=add_fit)
#merge df
com_add_df <- bind_cols(list(HK68_add,Bk79_add,Bei89_add,Mos99_add,Bris07_add,ND16_add)) %>%
  select(2:3,5,10,15,20,25,30)
textsize <- 7


p2  <- ggpairs(com_add_df,columns=3:8,lower =list(continuous=wrap(ggally_points,size=0.1,alpha=1)),
		upper = list(continuous="blank"),
		diag  = list(continuous ="blank")) +
  theme_cowplot(9) +
  theme(plot.title=element_text(size=textsize,face="bold"),
            panel.grid.major = element_blank(),
            legend.text=element_text(size=textsize,face="bold"),
            axis.title=element_text(size=textsize,face="bold"),
            axis.text=element_text(size=textsize,face="bold"))+
  theme(plot.title=element_text(size=textsize,face="bold"),
            panel.grid.major = element_blank(),
            legend.text=element_text(size=textsize,face="bold"),
            axis.title=element_text(size=textsize,face="bold"),
            axis.text=element_text(size=textsize,face="bold")) +
  scale_x_continuous(n.breaks = 3.5) +
  scale_y_continuous(n.breaks = 4)
ggsave('graph/HA_compare_ADD.png',p2,width = 4, height = 4)
cor_add <- cor(data.table(select(com_add_df, 3:8)))
print (cor(data.table(select(com_add_df, 3:8))))

#combine the epi cor and add cor
cor_combined <- cor_epi
cor_combined[upper.tri(cor_combined)] <- cor_add[upper.tri(cor_add)]

diag(cor_combined) <- NA # remove diagonal for better plotting (optional)
#reformat the data
StrainLevels <- c('HK68','Bk79','Bei89','Mos99','Bris07','ND16')
cor_df <- as.data.frame(cor_combined) %>%
	dplyr::mutate(strain1 = row.names(.)) %>%
  	tidyr::gather("strain2", "value", -strain1) %>%
	mutate(strain1=factor(strain1,levels=StrainLevels))%>%
	mutate(strain2=factor(strain2,levels=StrainLevels))
cor_df[,3]<-round(cor_df[,3], digits = 2)
#make heatmap plot
textsize    <- 7
heat_p <- ggplot(cor_df, aes(strain1, strain2)) +
	geom_tile(aes(fill = value), colour = "white",size =0.1) +
  	scale_fill_gradient2(low = "blue", high = "red", mid = "white",
						 midpoint = 0, limit = c(-1,1), space = "Lab",
						 name="Pearson\nCorrelation") +
  	scale_x_discrete(expand = c(0, 0)) +
  	scale_y_discrete(expand = c(0, 0)) +
  	labs(x=NULL, y=NULL) +
  	theme(legend.position = "none",
		  axis.ticks = element_blank(),
		  legend.text=element_text(size=textsize,face="bold"),
		  axis.text=element_text(size=textsize,face="bold",color='black'),
		  axis.text.x=element_text(angle=0,hjust=0.5,vjust=0.5),
		  axis.title=element_text(size=textsize,face="bold"))

#add correlation coefficients on the heatmap
heat_p <- heat_p+
	geom_text(aes(strain1, strain2, label = value), color = "black", size = 2) +
	theme(
		axis.title.x = element_blank(),
  		axis.title.y = element_blank(),
  		panel.grid.major = element_blank(),
  		panel.border = element_blank(),
  		panel.background = element_blank(),
  		axis.ticks = element_blank())
ggsave('graph/HA_heatmap_epi_add_cor.png',bg='white',heat_p,width = 2.3, height = 1.9)