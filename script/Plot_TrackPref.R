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

strain_levels <- c('HK68','Bk79','Bei89','Mos99','Vic11','HK19')
YearPrefTable <- read_tsv('result/Fits_ByYear.tsv') %>%
                   mutate(Background=factor(Background,levels=strain_levels))

palette  <- c(brewer.pal(9,"Set3"))
palette  <- c(palette[3:12])
textsize <- 7
p <-  ggplot(YearPrefTable, aes(x=Year,y=Meanfit,group=Background,color=Background)) +
        geom_line(size=0.4) + 
        geom_point(size=0.4) +
        geom_errorbar(aes(ymin=Meanfit-Stdfit, ymax=Meanfit+Stdfit), width=0) +
        scale_color_manual(values=palette) + 
        theme_cowplot(12) +
        theme(plot.title=element_text(size=textsize,face="bold"),
              axis.title=element_text(size=textsize,face="bold"),
              axis.text=element_text(size=textsize,face="bold"),
              legend.title=element_blank(),
              legend.text=element_text(size=textsize,face="bold"),
              legend.position='right') + 
        ylab(bquote(bold(Mean~Fitness)))
ggsave('graph/FitsByYear2.png',p,height=1.8,width=4.5)
