# Title     : plot cor of peak_epi and model_epi
# Objective : 
# Created by: yiquan
# Created on: 5/25/21
#R code
library(ggplot2)
library(scales)
library(ggpubr)
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


textsize <-7

peak_epi  <- read_csv('result/peak_epi.csv')
peak_epi <- peak_epi %>% filter(abs(peak_epi) > 0.01)%>%
  filter(!grepl('N|S|L|F|G|T', epi_resi))


#change the l1 to 328 and add mut1 and mut2 columnn
CompareEPI <- function (model_epi_path,peak_epi){
	model_epi  <- read_csv(model_epi_path)
	model_epi <- model_epi %>%
		mutate(L1=recode(L1, `1` = '328',`2` = '329', `3` = '344',`4` = '367',`5` = '368',`6` = '369',`7` = '370'),
			   L2=recode(L2, `1` = '328',`2` = '329', `3` = '344',`4` = '367',`5` = '368',`6` = '369',`7` = '370'))
	model_epi <- model_epi %>%
		mutate(mut1=paste(model_epi$AA1, model_epi$L1, sep = ''),
			   mut2=paste(model_epi$AA2, model_epi$L2, sep = ''))
	#add charge pair column
	model_epi <- model_epi %>%
		mutate(charge1=recode(AA1,'K'="+chg",'R'="+chg",'D'="-chg",'E'="-chg",.default ="neutral"),
			   charge2=recode(AA2,'K'="+chg",'R'="+chg",'D'="-chg",'E'="-chg",.default ="neutral"))
	model_epi <- model_epi %>%
		mutate(charge=paste(model_epi$charge1, model_epi$charge2, sep = '/')) %>%
		mutate(epi_resi=paste(model_epi$mut1, model_epi$mut2, sep = '-'))
	model_epi <- model_epi %>%
		mutate(charge=recode(model_epi$charge,"-chg/+chg"='+chg/-chg',"+chg/neutral"='neutral/+chg',"-chg/neutral"='neutral/-chg')) %>%
		filter(epi_resi %in% peak_epi$epi_resi) %>%
		select(epi_resi,EPI)

	#overlapping natural and model epistasis
	common_df <- merge(model_epi,peak_epi) %>%
		separate(epi_resi,c('aa1','aa2'),'-') %>%
		mutate(charge1=ifelse(grepl('K|R',aa1),'+',ifelse(grepl('D|E',aa1),'-','neutral')),
			   charge2=ifelse(grepl('K|R',aa2),'+',ifelse(grepl('D|E',aa2),'-','neutral')))%>%
		mutate(paircharge=ifelse(charge1 == 'neutral'|charge2=='neutral', 'Non_charge',ifelse((charge1 == '+'&charge2=='+')|(charge1 == '-'&charge2=='-'),'Same_charge','diff_charge'))) %>%
		select(aa1,EPI,aa2,peak_epi,paircharge)
	return(common_df)
}
PlotCorEPI <- function (common_df,title){
	p1 <- ggplot(common_df,aes(x=EPI,y=common_df$peak_epi,group=paircharge)) +
		geom_point(aes(shape=paircharge, color=paircharge),size=1) +
		theme_cowplot(12) +
		scale_color_manual(values=c("#E69F00", "#56B4E9")) +
		theme(plot.title=element_text(size=textsize,face="bold",hjust = 0.5),
			  axis.title=element_text(size=textsize,face="bold"),
			  axis.text=element_text(size=textsize,face="bold"),
			  legend.title=element_blank(),
			  legend.text=element_text(size=textsize,face="bold"),
			  legend.position='none') +
		ggtitle(title)+
		xlab(bquote(bold('Model EPI'))) +
		ylab(bquote(bold('Peak EPI')))
	cor(exp(common_df$EPI),exp(common_df$peak_epi))
	return(p1)
}



common_df1 <- CompareEPI('result/HK68_epi.csv',peak_epi)
common_df2 <- CompareEPI('result/Bk79_epi.csv',peak_epi)
common_df3 <- CompareEPI('result/Bei89_epi.csv',peak_epi)
common_df4 <- CompareEPI('result/Mos99_epi.csv',peak_epi)
common_df5 <- CompareEPI('result/Vic11_epi.csv',peak_epi)
common_df6 <- CompareEPI('result/HK19_epi.csv',peak_epi)

p1 <- PlotCorEPI(common_df1,'HK68')
p2 <- PlotCorEPI(common_df2,'Bk79')
p3 <- PlotCorEPI(common_df3,'Bei89')
p4 <- PlotCorEPI(common_df4,'Mos99')
p5 <- PlotCorEPI(common_df5,'Vic11')
p6 <- PlotCorEPI(common_df6,'HK19')
p <- ggarrange(p1,p2,p3,p4,p5,p6,nrow=3,ncol=2,common.legend = TRUE,legend="right")
ggsave('graph/compare_natural_epi.png',p,height=4,width=3.6)