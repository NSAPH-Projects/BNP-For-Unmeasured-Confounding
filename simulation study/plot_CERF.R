###################################################################################
# ---                   Plots results                ---
# ---                 SIMULATION STUDY               ----
###################################################################################

# library 
library(ggplot2)
library(GGally)
library(hrbrthemes)
library(mvtnorm)
library(invgamma)
library(BNPmix)
library(truncnorm)
library(ggrepel)
library(TeachingDemos)
library(fanplot)
library(ggfan)

# load data
load("results_model.RData")

###################################################################################
# ---     VARIABILITY BOUND      ----
###################################################################################

variability_bound<- function(data_sim,medians_x,medians_smooth,funct_int,title,num){
  
  points_x_i=medians_x
  points_y=funct_int(points_x_i,data_sim[[1]]$U)
  
  DB=as.data.frame(cbind(Treatment=rep(points_x_i, each=sample),
                         Outcome=c(medians_smooth)))
  
  pdf(paste0("plot_",title,".pdf"),width=6, height=5)
  gg<-ggplot(DB, aes(x=Treatment, y= Outcome)) + 
    geom_fan() + 
    geom_line(aes(rep(points_x_i,sample), rep(apply(medians_smooth,2,median),sample)),
              size=1, col="#0C2AE8") +
    geom_line(aes(rep(points_x_i,sample), rep(points_y,sample)),
              size=0.7, col="red") +
    theme_minimal()+ 
    scale_fill_gradient2(low="#0C2AE8", high="#DCF3F5", mid="#59C7FF", 
                         midpoint=0.5) +
    annotate("text", x=points_x_i[20], y=max(points_y,medians_smooth), 
             label=paste0("Scenario ",num,""), size=4.5)
  print(gg)
  dev.off()
}


# print plots
variability_bound(data_sim=data_sim_1,
                         medians_x=median_1s_4q$x_points,
                         medians_smooth=median_1s_4q$smooth02,
                         funct_int=fun_Y_1_int,
                         title="c1_4q", num=1)
variability_bound(data_sim=data_sim_2,
                         medians_x=median_2s_4q$x_points,
                         medians_smooth=median_2s_4q$smooth02,
                         funct_int=fun_Y_2_int,
                         title="c2_4q", num=2)
variability_bound(data_sim=data_sim_3,
                         medians_x=median_3s_4q$x_points,
                         medians_smooth=median_3s_4q$smooth02,
                         funct_int=fun_Y_3_int,
                         title="c3_4q", num=3)
variability_bound(data_sim=data_sim_4,
                         medians_x=median_4s_4q$x_points,
                         medians_smooth=median_4s_4q$smooth02,
                         funct_int=fun_Y_4_int,
                         title="c4_4q", num=4)


 # plots in appendix
variability_bound(data_sim=data_sim_1,
                  medians_x=median_1s_6q$x_points,
                  medians_smooth=median_1s_6q$smooth02,
                  funct_int=fun_Y_1_int,
                  title="c1_6q", num=1)
variability_bound(data_sim=data_sim_2,
                  medians_x=median_2s_6q$x_points,
                  medians_smooth=median_2s_6q$smooth02,
                  funct_int=fun_Y_2_int,
                  title="c2_6q", num=2)
variability_bound(data_sim=data_sim_3,
                  medians_x=median_3s_6q$x_points,
                  medians_smooth=median_3s_6q$smooth02,
                  funct_int=fun_Y_3_int,
                  title="c3_6q", num=3)
variability_bound(data_sim=data_sim_4,
                         medians_x=median_4s_6q$x_points,
                         medians_smooth=median_4s_6q$smooth02,
                         funct_int=fun_Y_4_int,
                         title="c4_6q", num=4)

###################################################################################
# ---     PLOTS SIMULATION STUDIES - simulation plots      ----
###################################################################################

funct_plot<-function(data_sim,funct_YU,funct_int,c){
  
  coeff=lm(data_sim$U ~ data_sim$X)$coeff
  
  pdf(paste0("sim_",c,".pdf"),width=6, height=6)
  gg<- ggplot(data_sim, aes(x=X, y=Y))+
    geom_function(fun = function(x) funct_int(x,data_sim$U), 
                  colour="red", lwd=1.5) +
    labs(y = "Outcome",
         x = "Treatment") +
    theme_minimal() +
    annotate("text", x=min(data_sim$X)+0.2, y=max(data_sim$Y)-0.7, 
             label=paste0("Scenario ",c,""), size=6)
  print(gg)
  dev.off()
}

funct_plot(data_sim=data_sim_1[[1]],
           funct_YU=fun_Y_1,
           funct_int=fun_Y_1_int,
           c=1)
funct_plot(data_sim=data_sim_2[[1]],
           funct_YU=fun_Y_2,
           funct_int=fun_Y_2_int,
           c=2)
funct_plot(data_sim=data_sim_3[[1]],
           funct_YU=fun_Y_3,
           funct_int=fun_Y_3_int,
           c=3)
funct_plot(data_sim=data_sim_4[[1]],
           funct_YU=fun_Y_4,
           funct_int=fun_Y_4_int,
           c=4)
