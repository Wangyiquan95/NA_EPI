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


PlotCompareFit_Rep <- function(table_data,strainname){
  textsize=7
  colorscale <- c(brewer.pal(9,"Set1"))
  table_data <- filter(table_data,strain==strainname)
  print (strainname)
  R1fit = log10(table_data$rep1_fit)
  R2fit = log10(table_data$rep2_fit)
  #R1fit[!is.finite(R1fit)] <- NA
  #R2fit[!is.finite(R2fit)] <- NA
  print (paste("Pearson Cor:", cor(R1fit,R2fit,use="complete.obs"),sep=' '))
  p <- ggplot(table_data,aes(x=log10(rep1_fit),y=log10(rep2_fit),color=mut_type)) +
	 geom_point(size=0.3) +
         theme_cowplot(12) +
         scale_color_manual(values=c('gray', colorscale)) +
	 theme(plot.title=element_text(size=textsize,face="bold",hjust = 0.5),
	       axis.title=element_text(size=textsize,face="bold"),
	       axis.text=element_text(size=textsize,face="bold"),
	       legend.title=element_blank(),
	       legend.text=element_text(size=textsize,face="bold"),
	       legend.position='none') +
	 xlab(bquote(bold('Replicate 1'))) +
	 ylab(bquote(bold('Replicate 2'))) +
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
p1 <- PlotCompareFit_Rep(table_data,'HK68')
p2 <- PlotCompareFit_Rep(table_data,'Bk79')
p3 <- PlotCompareFit_Rep(table_data,'Bei89')
p4 <- PlotCompareFit_Rep(table_data,'Mos99')
p5 <- PlotCompareFit_Rep(table_data,'Vic11')
p6 <- PlotCompareFit_Rep(table_data,'HK19')
p <- grid.arrange(p1,p2,p3,p4,p5,p6,nrow=3,ncol=2)
ggsave('graph/Compare_Rep.png',p,height=6,width=4)
