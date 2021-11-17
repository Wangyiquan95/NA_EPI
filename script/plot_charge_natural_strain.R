#R code
library(ggplot2)
library(scales)
library(RColorBrewer)
library(readr)
library(tidyr)
library(reshape)
library(stringr)
library(dplyr)
library(gridExtra)
require(cowplot)

plot_charge <- function(data_table, graphname){
  colorscale  <- c(brewer.pal(9,"Set1"))
  textsize <- 7
  p <- ggplot(data=data_table,aes(x=year,y=local_chg, color=local_chg)) +
         geom_point(size=0.1, alpha=0.1, shape=20,position='jitter') +
         #geom_ribbon(aes(ymax=local_chg+SD, ymin=local_chg-SD),colour = NA) +
         #scale_color_manual(values=colorscale,drop=FALSE) +
         #scale_fill_manual(values=alpha(colorscale,0.5),drop=FALSE) +

         theme_cowplot(12) +scale_color_gradientn(colours=c("blue", "gray90", "red"))+
         theme(legend.title=element_blank(),
               legend.key.size=unit(3,"mm"),
               legend.position='none',
               legend.justification='center',
               legend.text=element_text(size=textsize,face="bold"), 
               axis.text=element_text(size=textsize,face="bold"),
               axis.text.x=element_text(angle=90,hjust=1,vjust=0.5),
               axis.title=element_text(size=textsize,face="bold")) +
         ylab("Net charge") +
         xlab("Year") +
         scale_x_continuous(breaks=c(1968,1980,1990,2000,2010,2020),labels=c(1968,1980,1990,2000,2010,2020)) +
         scale_y_continuous(breaks=c(-4,-3,-2,-1,0,1,2,3,4),labels=c(-4,-3,-2,-1,0,1,2,3,4), limit=c(-4,4))
  ggsave(graphname,p,bg='white',height=2,width=2)
  }

data_NA <- read_tsv('result/HumanH3N2_NA_charge.tsv')
chglevels <- c(-3,-2,-1,0,1,2,3)

NA_table <- data.frame(table(select(data_NA, year,local_chg))) %>%
	group_by(year)%>%
	mutate(local_chg=factor(local_chg,levels=chglevels)) %>%
	mutate(freq_given_year = Freq / sum(Freq))

NA_table[is.na(NA_table)] <-0
NA_table$year <- as.numeric(as.character(NA_table$year))

p1 <- ggplot(NA_table,aes(x=year,y=freq_given_year,group=local_chg,color=local_chg)) +
	geom_line() +
	theme_cowplot(12) +scale_color_brewer(palette = "Set2")+
	theme(plot.title=element_text(size=7,face="bold",hjust = 0.5),
		  axis.title=element_text(size=7,face="bold"),
		  axis.text=element_text(size=7,face="bold",color = 'black'),
		  axis.text.x=element_text(angle=90,hjust=1,vjust=0.5),
		  legend.title=element_blank(),
		  legend.text=element_text(size=7,face="bold"),
		  legend.position='top',legend.key.size = unit(0.3, 'cm')) +
	xlab(bquote(bold('Year'))) +
	ylab(bquote(bold('')))+
	guides(color=guide_legend(nrow=2, byrow=TRUE))+
    scale_x_continuous(breaks=c(1968,1980,1990,2000,2010,2020),labels=c(1968,1980,1990,2000,2010,2020))+
	scale_y_continuous(limits=c(0,1),breaks=c(0,0.5,1),labels=c('0%','50%','100%'))

ggsave('graph/natural_charge_freq_by_year.png',bg='white',p1,width = 2, height = 2)
# plot_charge(data_NA, 'graph/NA_year_vs_charge.png')
