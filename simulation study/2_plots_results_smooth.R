###################################################################################
# ---                   Plots results                ---
# ---   SIMULATION PAPER: 300 Samples with 5000 ss   ----
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

#load results
setwd("C:/Users/dafne/Desktop/DOTTORATO/kate francesca/SIMULATIONS")
load("results_model_stand.RData")

###################################################################################
# ---     VARIABILITY BOUND      ----
###################################################################################

variability_bound<- function(data_sim,medians,funct_int,title){
  
  min_x=mean(unlist(sapply(1:sample, function(i) min(data_sim[[i]]$X))))
  max_x=mean(unlist(sapply(1:sample, function(i) max(data_sim[[i]]$X))))
  points_x_i=seq(min_x,max_x,length.out=200)
  points_y=funct_int(points_x_i,data_sim[[1]]$U)
  
  DB=as.data.frame(cbind(x_lab=rep(points_x_i, each=sample),
                         outcome=c(medians$with)))
  
  #setwd("C:/Users/dafne/Desktop/plot temporanei/UC-samples")
  #pdf(paste0("var_",title,"_intercept.pdf"),width=9, height=7)
  gg<-ggplot(DB, aes(x=x_lab, y= outcome)) + 
    geom_fan() + 
    #geom_interval(intervals=c(0))+
    #lines(points_x_i, apply(medians$with,2,median), col="red", type = "l", lwd = 4) +
    geom_line(aes(rep(points_x_i,sample), rep(apply(medians$with,2,median),sample)),
              size=1.1, col="#0C2AE8") +
    geom_line(aes(rep(points_x_i,sample), rep(points_y,sample)),
              size=1, col="red") +
    #scale_linetype_manual(values=c("solid"))+
    theme(
      panel.background = element_rect(fill='transparent'),
      plot.background = element_rect(fill ="#F8F1EE"),
      panel.grid.major = element_blank(),
      #panel.grid.minor = element_blank(),
      legend.background = element_rect(fill='transparent'),
      legend.box.background = element_rect(fill='transparent')
      #legend.position="bottom"
    )+ 
    scale_fill_gradient2(low="#0C2AE8", high="#DCF3F5", mid="#59C7FF", 
                         midpoint=0.5) +
    scale_color_manual(name = " ", 
                       values = c("median" = "#0C2AE8", "real" = "#E69F00"))
  print(gg)
  #dev.off()
}

variability_bound(data_sim=data_sim_1,
                  medians=median_1,
                  funct_int=fun_Y_1_int,
                  title="c1_")

variability_bound(data_sim=data_sim_2,
                  medians=median_2,
                  funct_int=fun_Y_2_int,
                  title="c2_")


variability_bound(data_sim=data_sim_3,
                  medians=median_3s_4q,
                  funct_int=fun_Y_3_int,
                  title="c3_4q")
variability_bound(data_sim=data_sim_3,
                  medians=median_3s_6q,
                  funct_int=fun_Y_3_int,
                  title="c3_6q")

variability_bound(data_sim=data_sim_4,
                  medians=median_4s_4q,
                  funct_int=fun_Y_4_int,
                  title="c4_4q")
variability_bound(data_sim=data_sim_4,
                  medians=median_4s_6q,
                  funct_int=fun_Y_4_int,
                  title="c4_6q")

###################################################################################
# ---     VARIABILITY BOUND -- white background for paper      ----
###################################################################################

variability_bound<- function(data_sim,medians,funct_int,title,num){
  
  min_x=mean(unlist(sapply(1:sample, function(i) min(data_sim[[i]]$X))))
  max_x=mean(unlist(sapply(1:sample, function(i) max(data_sim[[i]]$X))))
  points_x_i=seq(min_x,max_x,length.out=200)
  points_y=funct_int(points_x_i,data_sim[[1]]$U)
  
  DB=as.data.frame(cbind(Treatment=rep(points_x_i, each=sample),
                         Outcome=c(medians$with)))
  
  setwd("C:/Users/dafne/Dropbox/unmeasure confounder/CODE/SIMULATION section - paper")
  pdf(paste0("plot_",title,"_intercept.pdf"),width=9, height=7)
  gg<-ggplot(DB, aes(x=Treatment, y= Outcome)) + 
    geom_fan() + 
    #geom_interval(intervals=c(0))+
    #lines(points_x_i, apply(medians$with,2,median), col="red", type = "l", lwd = 4) +
    geom_line(aes(rep(points_x_i,sample), rep(apply(medians$with,2,median),sample)),
              size=1.1, col="#0C2AE8") +
    geom_line(aes(rep(points_x_i,sample), rep(points_y,sample)),
              size=1, col="red") +
    #scale_linetype_manual(values=c("solid"))+
    theme_minimal()+ 
    scale_fill_gradient2(low="#0C2AE8", high="#DCF3F5", mid="#59C7FF", 
                         midpoint=0.5) +
    scale_color_manual(name = " ", 
                       values = c("median" = "#0C2AE8", "real" = "#E69F00"))+
    annotate("text", x=points_x_i[10], y=max(points_y,medians$with), 
             label=paste0("Scenario ",num,""), size=6)
  print(gg)
  dev.off()
}

variability_bound(data_sim=data_sim_1,
                  medians=median_1s_4q,
                  funct_int=fun_Y_1_int,
                  title="c1_4q",
                  num=1)
variability_bound(data_sim=data_sim_1,
                  medians=median_1s_6q,
                  funct_int=fun_Y_1_int,
                  title="c1_6q",
                  num=1)

variability_bound(data_sim=data_sim_2,
                  medians=median_2s_4q,
                  funct_int=fun_Y_2_int,
                  title="c2_4q",
                  num=2)
variability_bound(data_sim=data_sim_2,
                  medians=median_2s_6q,
                  funct_int=fun_Y_2_int,
                  title="c2_6q",
                  num=2)

variability_bound(data_sim=data_sim_3,
                  medians=median_3s_4q,
                  funct_int=fun_Y_3_int,
                  title="c3_4q",
                  num=3)
variability_bound(data_sim=data_sim_3,
                  medians=median_3s_6q,
                  funct_int=fun_Y_3_int,
                  title="c3_6q",
                  num=3)

variability_bound(data_sim=data_sim_4,
                  medians=median_4s_4q,
                  funct_int=fun_Y_4_int,
                  title="c4_4q",
                  num=4)
variability_bound(data_sim=data_sim_4,
                  medians=median_4s_6q,
                  funct_int=fun_Y_4_int,
                  title="c4_6q",
                  num=4)

###################################################################################
# ---     PLOTS SIMULATION STUDIES - white background for paper      ----
###################################################################################

funct_plot<-function(data_sim,funct_YU,funct_int,c){
  
  coeff=lm(data_sim$U ~ data_sim$X)$coeff
  
  setwd("C:/Users/dafne/Dropbox/unmeasure confounder/CODE/SIMULATION section - paper")
  pdf(paste0("sim_",c,".pdf"),width=6, height=6)
  gg<- ggplot(data_sim, aes(x=X, y=Y))+
    #geom_function(fun = function(x) funct_YU(x,coeff[1]+coeff[2]*x), 
    #              colour="black", lwd=1.5, lty=4) +
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
