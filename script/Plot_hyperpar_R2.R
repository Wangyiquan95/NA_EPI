# Title     : plot reg/lr/bs affect model prediction R
# Objective :
# Created by: yiquan
# Created on: 6/15/21

library(readxl)
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
library(ggpubr)


StrainLevels <- c('HK68','Bk79','Bei89','Mos99','Vic11','HK19')
colorscale <- c(brewer.pal(8,"Set2"))
textsize <- 7
plot_corr <- function(df,x,y,strain,xname,yname){
  p <- ggplot(df,aes(x,y,group=strain,color = strain))+
    theme_cowplot(9) +geom_smooth(size=0.5,se = FALSE)+
    scale_color_manual(values=colorscale) + ylim(-0.3,1.3)+
    xlab(xname) + ylab(yname)+ labs(color='')+
    theme(plot.title=element_text(size=textsize,face="bold"),
          panel.grid.major = element_blank(),
          legend.position = "top",
          legend.text=element_text(size=textsize,face="bold"),
          legend.key.size = unit(0.5,"line"),
          axis.text=element_text(size=textsize,face="bold"),
          axis.title=element_text(size=textsize,face="bold",color='black'))
  return(p)
}
df <- read_excel('result/hyperparameter_R2.xlsx')
lr_df <- df %>%
  select(strain1,lr,add1,epi1) %>%
  group_by(strain1)  %>% drop_na()%>%
  mutate(strain=factor(strain1,levels=StrainLevels))
reg_df <- df %>%
  select(strain2,reg,add2,epi2) %>%
  group_by(strain2)%>% drop_na() %>%
  mutate(strain=factor(strain2,levels=StrainLevels))
bs_df <- df %>%
  select(strain3,bs,add3,epi3) %>%
  group_by(strain3)%>% drop_na() %>%
  mutate(strain=factor(strain3,levels=StrainLevels))

lr_add <- plot_corr(lr_df,lr_df$lr,lr_df$add1,lr_df$strain,'Learning rate',"Pearson Correlation")
reg_add <- plot_corr(reg_df,reg_df$reg,reg_df$add2,reg_df$strain,'Regularization',"Pearson Correlation")
bs_add <- plot_corr(bs_df,bs_df$bs,bs_df$add3,bs_df$strain,'Batch size',"Pearson Correlation")

lr_epi <- plot_corr(lr_df,lr_df$lr,lr_df$epi1,lr_df$strain,'Learning rate',"Pearson Correlations")
reg_epi <- plot_corr(reg_df,reg_df$reg,reg_df$epi2,reg_df$strain,'Regularization',"Pearson Correlation")
bs_epi <- plot_corr(bs_df,bs_df$bs,bs_df$epi3,bs_df$strain,'Batch size',"Pearson Correlation")

add <- ggarrange(lr_add,lr_epi,reg_add,reg_epi,bs_add,bs_epi,nrow=3,ncol=2,common.legend = TRUE,legend="right")
ggsave('graph/hyperpar_r.png',bg='white',add,height=4.0,width=5)

###plot regularization vs prediction R2
reg_r2_df <- read_csv('result/reg_r2.csv')
add_model_reg_r2_df <-read_csv('result/add_only_model_reg_R2.csv')%>%
  group_by(strain)%>%
  mutate(strain=factor(strain,levels=StrainLevels))

add_reg_r2 <- plot_corr(add_model_reg_r2_df,add_model_reg_r2_df$reg,add_model_reg_r2_df$R2_mean,add_model_reg_r2_df$strain,'Regularization',expression(bold('R'^2)))+
  scale_x_log10()+geom_errorbar(aes(ymin=add_model_reg_r2_df$R2_mean-add_model_reg_r2_df$R2_std, ymax=add_model_reg_r2_df$R2_mean+add_model_reg_r2_df$R2_std), width=0)
ggsave('graph/add_reg_r2.png',bg='white',add_reg_r2,height=2.0,width=5)
reg_r2_df <- reg_r2_df %>%
  group_by(strain)%>%
  mutate(strain=factor(strain,levels=StrainLevels))


reg_r2 <- plot_corr(reg_r2_df,reg_r2_df$reg,reg_r2_df$R2_mean,reg_r2_df$strain,'Regularization',expression(bold('R'^2)))+
  scale_x_log10()+geom_errorbar(aes(ymin=reg_r2_df$R2_mean-reg_r2_df$R2_std, ymax=reg_r2_df$R2_mean+reg_r2_df$R2_std), width=0)

reg_var <- plot_corr(reg_r2_df,reg_r2_df$reg,reg_r2_df$Var_mean,reg_r2_df$strain,'Regularization','Variance')+
  scale_x_log10()+geom_errorbar(aes(ymin=reg_r2_df$Var_mean-reg_r2_df$Var_std, ymax=reg_r2_df$Var_mean+reg_r2_df$Var_std), width=0)+
  ylim(-0.1,2)

reg_ran_var <- plot_corr(reg_r2_df,reg_r2_df$reg,reg_r2_df$Ran_Var_mean,reg_r2_df$strain,'Regularization','Variance(S2)')+
  scale_x_log10()+geom_errorbar(aes(ymin=reg_r2_df$Ran_Var_mean-reg_r2_df$Ran_Var_std, ymax=reg_r2_df$Ran_Var_mean+reg_r2_df$Ran_Var_std), width=0)+
  ylim(0,1)
r2 <- ggarrange(reg_r2,reg_var,nrow=2,ncol=1,common.legend = TRUE,legend="right")
ggsave('graph/reg_r2.png',bg='white',r2,height=4.0,width=5)
