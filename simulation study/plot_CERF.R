###################################################################################
# ---                 Plots Functions                ---
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

###################################################################################
# ---     VARIABILITY BOUND      ----
###################################################################################

ggplot_posterior <- function(data_sim,
                             model_adj,
                             x_split_method = "split_4quantile_x",
                             probs = NULL,
                             pts = NULL,
                             n_group = 10,
                             funct_YU,
                             bandwidth = 0.2,
                             smooth_curve = TRUE,
                             num) {
  
  
  # causal effect for a grid of points
  points_x_i = seq(min(data_sim$X), max(data_sim$X), length.out = 200)
  points_y_iteractions = sapply(points_x_i, function(x)
    curve_ADJ(
      x,
      post_chain = model_adj,
      training_data = data_sim,
      x_split_method =  x_split_method,
      probs = probs,
      pts = pts,
      n_group = n_group
    ))
  
  # plot
  par(mfrow = c(1, 1), bg = "white")
  # posterior median
  if (smooth_curve) {
    
    smooth_values <- sapply(1:dim(points_y_iteractions)[1], function(i) {
      ksmooth(points_x_i, points_y_iteractions[i, ],  "normal", bandwidth = bandwidth)$y
    })
    y_est = apply(t(smooth_values), 2, median, na.rm = TRUE)
    lower95 <- apply(t(smooth_values), 2, quantile, prob = 0.05)
    upper95 <- apply(t(smooth_values), 2, quantile, prob = 0.95)
    
  } else{
    
    y_est <- apply(points_y_iteractions, 2, median, na.rm = TRUE)
    lower95 <- apply(points_y_iteractions, 2, quantile, prob = 0.05)
    upper95 <- apply(points_y_iteractions, 2, quantile, prob = 0.95)
    
  }
  
  points_y_i <- funct_YU(points_x_i, data_sim$U)
  DB = as.data.frame(cbind(Treatment = data_sim$X,
                           Outcome = data_sim$Y))
  DB2 = as.data.frame(
    cbind(
      Treatment = points_x_i,
      Outcome = y_est,
      points_y_i = points_y_i,
      lower95 = lower95,
      upper95 = upper95
    )
  )
  gg <- ggplot() +
    geom_point(
      data = DB,
      aes(x = Treatment, y = Outcome),
      size = 0.3,
      col = "#D3D3D3"
    ) +
    geom_line(
      data = DB2,
      aes(x = Treatment, y = Outcome),
      size = 1.5,
      col = "#0C2AE8"
    ) +
    geom_line(
      data = DB2,
      aes(Treatment, points_y_i),
      size = 1,
      col = "red"
    ) +
    geom_ribbon(
      data = DB2,
      aes(x = Treatment, ymin = lower95, ymax = upper95),
      color = "#59C7FF",
      fill = "#59C7FF",
      alpha = 0.5
    ) +
    theme(
      panel.background = element_rect(fill = 'white'),
      plot.background = element_rect(fill = "white"),
      axis.title = element_text(size = 14),
      legend.text = element_text(size = 10),
      legend.title = element_text(size = 12),
      plot.title = element_text(hjust = 0.5),
      axis.text.x = element_text(size = 14),
      axis.text.y = element_text(size = 14),
      title = element_text(size = 14),
      legend.background = element_rect(fill = 'transparent'),
      panel.grid.major = element_line(color = "grey", size = 0.1)
    ) +
    ggtitle(paste0("Scenario ", num, ""))
  print(gg)
}


plot_posterior<- function(data_sim,post_chian,funct_YU){
  
  
  # causal effect for a grid of points
  points_x_i=seq(min(data_sim$X),max(data_sim$X),length.out=200)
  points_y_iteractions=sapply(points_x_i, function(x) curve_chains(x,data_sim, post_chian))
  quantile_adj=apply( points_y_iteractions,2,quantile, prob=c(0.05,0.95))
  
  # plot
  par(mfrow=c(1,1),bg = "white")
  # raw data
  plot(data_sim$X,data_sim$Y,pch=19, cex=0.1, xlab="Treatment", ylab="Outcome", cex.lab =1.5, cex.axis =1.5, axes=T,  col = "#D3D3D3", 
  )
  # posterior median
  points(points_x_i,apply(points_y_iteractions, 2, median, na.rm=TRUE),  cex=0.5, pch=16, lwd = 0.5,col=rainbow(7)[6]) 
  # true
  curve(funct_YU(x,data_sim$U), min(data_sim$X),max(data_sim$X), col="#900000",lwd=2, add=TRUE)
  # 95% CI
  polygon(c(points_x_i, rev(points_x_i)), 
          c(quantile_adj[1,], rev(quantile_adj[2,])),
          col = "#0078EF7D",border = "#0078EF7D")
  grid(nx = NULL, ny = NULL,
       lty = 2,      # Grid line type
       col = "gray", # Grid line color
       lwd = 2)
}

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
