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
library(data.table)
library(GGally)
library(e1071) 
library(sinaplot)
library(ggforce)
require(cowplot)

plot_fit_sina <- function(table_fit, StrainLevels){ 
  textsize <- 7
  table_fit <- table_fit %>%
                  mutate(strain=factor(strain, levels=rev(StrainLevels)))
  colorscale <- c(brewer.pal(9,"Set1"))
  p <- ggplot(table_fit, aes(x=strain, y=fit, color=mut_type, group=strain)) +
         geom_sina(size=0.2) +
         theme_cowplot(12) +
         scale_color_manual(values=c('gray', colorscale)) +
         theme(plot.title=element_text(size=textsize,face="bold"),
               panel.grid.major = element_blank(),
               legend.position = "none",
               legend.text=element_text(size=textsize,face="bold"),
               legend.key.size = unit(0.5,"line"),
               axis.title=element_blank(),
               axis.text=element_text(size=textsize,face="bold")) +
         coord_flip()
  ggsave('graph/Lib_fit_sina.png',p,width=2,height=3.5)
  }

plot_pairs <- function(table_fit_cast){
  textsize <- 7
  colorscale <- c(brewer.pal(9,"Set1"))
  table_fit_cast <- table_fit_cast 
  p  <- ggpairs(data=table_fit_cast, columns=3:8,
                mapping=ggplot2::aes(colour=mut_type),
		lower =list(continuous=wrap(ggally_points,size=0.1,alpha=1)),
		upper = list(continuous="blank"),
		diag  = list(continuous ="blank")) +
          theme_cowplot(12) +
	  theme(plot.title=element_text(size=textsize,face="bold"),
            panel.grid.major = element_blank(),
            legend.text=element_text(size=textsize,face="bold"),
            axis.title=element_text(size=textsize,face="bold"),
            axis.text=element_text(size=textsize,face="bold")) +
      scale_x_continuous(n.breaks = 4) +
      scale_y_continuous(n.breaks = 4)
  for(i in 1:p$nrow) {
    for(j in 1:p$ncol){
      p[i,j] <- p[i,j] + 
        scale_fill_manual(values=c('gray', colorscale)) +
        scale_color_manual(values=c('gray', colorscale))
      }
    }
  ggsave('graph/LibCorPairs.png',p,width=4.4,height=4.4)
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

sizing <- function(mut_type){
  if (mut_type=='mut'){return (1)}
  else {return (3)}
  }

StrainLevels <- c('HK68','Bk79','Bei89','Mos99','Vic11','HK19')
table_fit <- read_tsv('result/NA_compile_results.tsv') %>%
                mutate(fit=log10(fit)) %>%
                mutate(strain=factor(strain,levels=StrainLevels)) %>%
                mutate(mut_type=mapply(coloring, ID)) %>%
                mutate(size=mapply(sizing, mut_type)) %>%
                mutate(mut_type=factor(mut_type, levels=c('mut', StrainLevels))) %>%
                arrange(mut_type)
table_fit_cast <- cast(table_fit,ID+mut_type~strain,value='fit') %>%
                     arrange(mut_type)
print (cor(data.table(select(table_fit_cast, HK68, Bk79, Bei89, Mos99, Vic11, HK19))))
plot_fit_sina(table_fit, StrainLevels)
plot_pairs(table_fit_cast)

#Calculate the mean fitness
fit_sum <-  table_fit %>%
  group_by(strain) %>% summarise(avg_fit=mean(fit))