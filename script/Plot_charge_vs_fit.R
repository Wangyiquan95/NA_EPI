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


plot_pI_vs_fit <- function(table_data,strainname){
  textsize=7
  colorscale <- c(brewer.pal(12,"Set3"))
  table_data <- filter(table_data,strain==strainname)
  print (strainname)
  p <- ggplot(table_data,aes(x=charge,y=log10(fit))) +
	 geom_point(size=0.3, color='gray70', alpha=0.1)+#,shape=20) + 
         geom_smooth(method="loess", se=FALSE, fullrange=FALSE, level=0.95, color=colorscale[1] ,size=0.6) +
         theme_cowplot(12) +
         #scale_fill_manual(values=c('gray', colorscale)) +
	 theme(plot.title=element_text(size=textsize,face="bold",hjust = 0.5),
	       axis.title=element_text(size=textsize,face="bold"),
	       axis.text=element_text(size=textsize,face="bold"),
	       legend.title=element_blank(),
	       legend.text=element_text(size=textsize,face="bold"),
	       legend.position='none') +
	 xlab(bquote(bold('Net charge'))) +
	 ylab('Fitness') +
         ggtitle(strainname)
  return (p)
  }

coloring <- function(ID){
  if (ID=='NDRSKDL'){return ('HK68')}
  else if (ID=='KNKSEES'){return ('Bk79')}
  else if (ID=='KNKGEEL'){return ('Bei89')}
  else if (ID=='KNESEKL'){return ('Mos99')}
  else if (ID=='KTENETS'){return ('Vic11')}
  else if (ID=='KSENETS'){return ('HK19')}
  else {return ('mut')}
  }

StrainLevels <- c('HK68','Bk79','Bei89','Mos99','Vic11','HK19')
table_data <- read_tsv('result/NA_compile_results.tsv') %>%
                mutate(mut_type=mapply(coloring, ID)) %>%
                mutate(mut_type=factor(mut_type, levels=c('mut', StrainLevels))) %>%
                arrange(mut_type)
p1 <- plot_pI_vs_fit(table_data,'HK68')
p2 <- plot_pI_vs_fit(table_data,'Bk79')
p3 <- plot_pI_vs_fit(table_data,'Bei89')
p4 <- plot_pI_vs_fit(table_data,'Mos99')
p5 <- plot_pI_vs_fit(table_data,'Vic11')
p6 <- plot_pI_vs_fit(table_data,'HK19')
p <- grid.arrange(p1,p2,p3,p4,p5,p6,nrow=3,ncol=2)
ggsave('graph/Compare_charge_vs_fit_bK79.png',bg='white',p2,height=2,width=2)
ggsave('graph/Compare_charge_vs_fit.png',bg='white',p,height=6,width=4)
