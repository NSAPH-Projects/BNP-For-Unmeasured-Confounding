############## PLOTS ###############################

load("~/chains_APPLICATION.RData")

####################################################


pdf("results_medan.pdf",width=9, height=7)
par(mfrow = c(1, 1))
plot(
  data_analysis$X,
  data_analysis$Y,
  pch = 19,
  cex = 0.3,
  col = 'dark grey',
  xlab = "PM2.5 exposure",
  ylab = "log hospitalization rate",
  ylim = c(-3.5,-1.8),
  xlim = c(min(estimated_CERF$points_x_i)+0.55,max(estimated_CERF$points_x_i)-0.55)
)
lines(estimated_CERF$points_x_i, estimated_CERF$smooth_YX, col = rainbow(7)[1], lwd = 2)
lines(estimated_CERF$points_x_i, estimated_CERF$smooth_YXU, col = rainbow(7)[3], lwd = 2)
lines(estimated_CERF$points_x_i, estimated_CERF$smooth_adj, col = rainbow(7)[6], lwd = 2)
legend(
  "topright",
  legend = c("YX", "YXU", "BNP-NC"),
  col = rainbow(7)[c(1, 3, 6)],
  lwd = 2
)
dev.off()

quantile_YX=apply(estimated_CERF$points_YX,2,quantile, prob=c(0.05,0.95))
quantile_YX=apply(estimated_CERF$points_YX,2,quantile, prob=c(0.025,0.975))
quantile_YXU=apply(estimated_CERF$points_YXU,2,quantile, prob=c(0.05,0.95))
quantile_adj=apply(estimated_CERF$points_adj,2,quantile, prob=c(0.05,0.95))

smooth_funct<-function(values){
  return(smooth_adj = ksmooth(
    estimated_CERF$points_x_i, values,
    "normal",bandwidth = 0.5)$y)
}

pdf("results_90CI.pdf",width=9, height=7)
par(mfrow = c(1, 1))
plot(
  data_analysis$X,
  data_analysis$Y,
  pch = 19,
  cex = 0.3,
  col = 'dark grey',
  xlab = "PM2.5 exposure",
  ylab = "log hospitalization rate",
  ylim = c(-3.5,-1.8),
  xlim = c(min(estimated_CERF$points_x_i)+0.55,max(estimated_CERF$points_x_i)-0.55)
)
polygon(c(estimated_CERF$points_x_i, rev(estimated_CERF$points_x_i)), 
        c(smooth_funct(quantile_YX[1,]), rev(smooth_funct(quantile_YX[2,]))),
        col = "#F8710C7D", border = "#F8710C7D")
polygon(c(estimated_CERF$points_x_i, rev(estimated_CERF$points_x_i)), 
        c(smooth_funct(quantile_YXU[1,]), rev(smooth_funct(quantile_YXU[2,]))),
        col = "#0BD97F7D",border = "#0BD97F7D")
polygon(c(estimated_CERF$points_x_i, rev(estimated_CERF$points_x_i)), 
        c(smooth_funct(quantile_adj[1,]), rev(smooth_funct(quantile_adj[2,]))),
        col = "#0078EF7D",border = "#0078EF7D")
lines(estimated_CERF$points_x_i, estimated_CERF$smooth_YX, col = rainbow(7)[1], lwd = 2)
lines(estimated_CERF$points_x_i, estimated_CERF$smooth_YXU, col = rainbow(7)[3], lwd = 2)
lines(estimated_CERF$points_x_i, estimated_CERF$smooth_adj, col = rainbow(7)[6], lwd = 2)
legend(
  "topright",
  legend = c("YX: median", "YXU: median", "BNP-NC: median",
             "YX: 90% CI", "YXU: 90% CI", "BNP-NC: 90% CI"),
  col = c(rainbow(7)[c(1, 3, 6)],"#F8710C7D","#0BD97F7D","#0078EF7D"),
  lwd = c(2,2,2,4,4,4)
)
dev.off()


quantile_YX=apply(estimated_CERF$points_YX,2,quantile, prob=c(0.025,0.975))
quantile_YXU=apply(estimated_CERF$points_YXU,2,quantile, prob=c(0.025,0.975))
quantile_adj=apply(estimated_CERF$points_adj,2,quantile, prob=c(0.025,0.975))


pdf("results_95CI.pdf",width=9, height=7)
par(mfrow = c(1, 1))
plot(
  data_analysis$X,
  data_analysis$Y,
  pch = 19,
  cex = 0.3,
  col = 'dark grey',
  xlab = "PM2.5 exposure",
  ylab = "log hospitalization rate",
  ylim = c(-3.5,-1.8),
  xlim = c(min(estimated_CERF$points_x_i)+0.55,max(estimated_CERF$points_x_i)-0.55)
)
polygon(c(estimated_CERF$points_x_i, rev(estimated_CERF$points_x_i)), 
        c(smooth_funct(quantile_YX[1,]), rev(smooth_funct(quantile_YX[2,]))),
        col = "#F8710C7D", border = "#F8710C7D")
polygon(c(estimated_CERF$points_x_i, rev(estimated_CERF$points_x_i)), 
        c(smooth_funct(quantile_YXU[1,]), rev(smooth_funct(quantile_YXU[2,]))),
        col = "#0BD97F7D",border = "#0BD97F7D")
polygon(c(estimated_CERF$points_x_i, rev(estimated_CERF$points_x_i)), 
        c(smooth_funct(quantile_adj[1,]), rev(smooth_funct(quantile_adj[2,]))),
        col = "#0078EF7D",border = "#0078EF7D")
lines(estimated_CERF$points_x_i, estimated_CERF$smooth_YX, col = rainbow(7)[1], 
      lwd = 2, lty = 2)
lines(estimated_CERF$points_x_i, estimated_CERF$smooth_YXU, col = rainbow(7)[3], 
      lwd = 2, lty = 2)
lines(estimated_CERF$points_x_i, estimated_CERF$smooth_adj, col = rainbow(7)[6], 
      lwd = 2, lty = 2)
legend(
  "topright",
  legend = c("YX: median", "YXU: median", "BNP-NC: median",
             "YX: 95% CI", "YXU: 95% CI", "BNP-NC: 95% CI"),
  col = c(rainbow(7)[c(1, 3, 6)],"#F8710C7D","#0BD97F7D","#0078EF7D"),
  lwd = c(2,2,2,4,4,4),
  lty = c(2,2,2,1,1,1)
)
dev.off()

