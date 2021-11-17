# Title     : plot trajectory
# Objective :
# Created by: yiquan
# Created on: 10/10/21
library(ggplot2)
library(dplyr)
library(readr)
library(RColorBrewer)
require(cowplot)
library(gghighlight)
library(ggpubr)

int_breaks <- function(x, n = 5) {
  l <- pretty(x, n)
  l[abs(l %% 1) < .Machine$double.eps ^ 0.5]
}
plot_trajectory <- function(df,strain,x_lab,y_lab){
  p <- ggplot(df,aes(x=step,y=E_score,group=Trajectory,colour = Trajectory))+
    geom_line()+gghighlight(trajectory_type==strain,label_params = list(size = 1.8))+#use_direct_label=FALSE)+
    theme_cowplot(12) +
    theme(plot.title=element_text(size=7,face="bold"),
          panel.grid.major = element_blank(),
          legend.position = "none",
          legend.text=element_text(size=7,face="bold"),
          legend.key.size = unit(0.5,"line"),
          axis.text=element_text(size=7,face="bold"),
          axis.title=element_text(size=7,face="bold"))+labs(x=x_lab,y = y_lab)+scale_x_continuous(breaks = int_breaks)
  return(p)
}

coloring <- function(ID){
  if (ID=='NDRSKDL-NNRSKDL-NNRSKDS-NNRSEDS-NNKSEDS'){return ('HK68')}
  else if (ID=='KNKSEES-KNKGEEL'){return ('Bk79')}
  else if (ID=='KNKGEEL-KNKSEEL-KNESEKL'){return ('Bei89')}
  else if (ID=='KNESEKL-KNESEKS-KNENETS'){return ('Mos99')}
  else if (ID=='KNENETS-KSENETS'){return ('Vic11')}
  else if (ID=='KSENETS-KSKNETS'){return ('HK19')}
  else {return ('mut')}
  }

HK68_df <- read_tsv('result/trajectory_prediction_HK68-4steps.tsv')%>%
                mutate(trajectory_type=mapply(coloring, Trajectory))
HK68_p <- plot_trajectory(HK68_df,'HK68','HK68 Trajectory','Fitness')

Bk79_df <- read_tsv('result/trajectory_prediction_Bk79-2steps.tsv')%>%
                mutate(trajectory_type=mapply(coloring, Trajectory))
Bk79_p <- plot_trajectory(Bk79_df,'Bk79','Bk79 Trajectory','Fitness')

Bei89_df <- read_tsv('result/trajectory_prediction_Bei89-3steps.tsv')%>%
                mutate(trajectory_type=mapply(coloring, Trajectory))
Bei89_p <- plot_trajectory(Bei89_df,'Bei89','Bei89 Trajectory','Fitness')

Mos99_df <- read_tsv('result/trajectory_prediction_Mos99-3steps.tsv')%>%
                mutate(trajectory_type=mapply(coloring, Trajectory))
Mos99_p <- plot_trajectory(Mos99_df,'Mos99','Mos99 Trajectory','Fitness')

Vic11_df <- read_tsv('result/trajectory_prediction_Vic11-1step.tsv')%>%
                mutate(trajectory_type=mapply(coloring, Trajectory))
Vic11_p <- plot_trajectory(Vic11_df,'Vic11','Vic11 Trajectory','Fitness')
HK19_df <- read_tsv('result/trajectory_prediction_HK19-1step.tsv')%>%
                mutate(trajectory_type=mapply(coloring, Trajectory))
HK19_p <- plot_trajectory(HK19_df,'HK19','HK19 Trajectory','Fitness')

trajectory_plot <- ggarrange(HK68_p,Bk79_p,Bei89_p,Mos99_p,Vic11_p,HK19_p,nrow=3,ncol=2)
ggsave('graph/HK68_trajectory.png',HK68_p,height=2.5,width=2.5,bg='white')
ggsave('graph/trajectory.png',trajectory_plot,height=7.5,width=5,bg='white')