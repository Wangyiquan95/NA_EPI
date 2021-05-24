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
library(gridExtra)
require(cowplot)

plotaafreq <- function(aatable,poi,xlab,ylab){
  aalevels <- c('R','K','D','E','L','N','S','T','G')
  postable    <- aatable %>%
                   filter(grepl(poi,mut)) %>%
                   mutate(aa=factor(aa,levels=aalevels))
  colorscale  <- c(brewer.pal(12,"Set3"))
  palette     <- c(colorscale[1],colorscale[3:12])
  textsize    <- 9
  p <- ggplot(postable,aes(year,freq,group=aa,color=aa)) +
         geom_line() +
         scale_color_manual(values=palette,drop=FALSE) +
         theme_cowplot() +
         theme(legend.title=element_blank(),
               #legend.key.size=unit(3,"mm"),
               legend.text=element_text(size=textsize,face="bold"),
               axis.text=element_text(size=textsize,face="bold"),
               axis.text.x=element_text(angle=0,hjust=0.5,vjust=0.5),
               axis.title=element_text(size=textsize,face="bold")) +
         ylab(ylab) +
         xlab(xlab) +
         scale_y_continuous(breaks=c(0,0.5,1),labels=c('0%','50%','100%')) +
         scale_x_continuous(breaks=c(1968,1980,1990,2000,2010,2020),labels=c(1968,1980,1990,2000,2010,2020)) +
         guides(color=guide_legend(ncol=1))
  return (p)
  }

aatable <- read_tsv('result/HumanN2Sweep_All.tsv') %>%
             mutate(aa=mapply(function(s){return(str_sub(s,-1,-1))},mut)) %>%
             mutate(mut=mapply(function(s){return(paste(str_sub(s,-1,-1),str_sub(s,2,-2),sep=''))},mut))
p_328 <- plotaafreq(aatable,'328',NULL,'residue 328')
p_329 <- plotaafreq(aatable,'329',NULL,'residue 329')
p_344 <- plotaafreq(aatable,'344',NULL,'residue 344')
p_367 <- plotaafreq(aatable,'367',NULL,'residue 367')
p_368 <- plotaafreq(aatable,'368',NULL,'residue 368')
p_369 <- plotaafreq(aatable,'369',NULL,'residue 369')
p_370 <- plotaafreq(aatable,'370',NULL,'residue 370')

#p <- grid.arrange(p_190,p_145,p_227,nrow=3)
p <- ggarrange(p_328,p_329,p_344,p_367,p_368,p_369,p_370,nrow=7,ncol=1,common.legend = TRUE,legend="right")
ggsave('graph/NatMutFreq_roi.png',p,height=6.4,width=5)
