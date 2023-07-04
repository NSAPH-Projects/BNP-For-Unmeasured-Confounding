###################################################################################
# ---          DATA ANALYSIS FOR THE PAPER         ---
# ---                    results                   ---
###################################################################################

# load results from BNP model
load("/n/dominici_nsaph_l3/projects/pm25-cardiovasculardisease-bnp/results_models_data_2.RData")

###################################################################################

plot(data_analysis$X,data_analysis$Y, pch=19, cex=0.5,
     xlab="treatment: X", ylab="outcome: Y", ylim=c(-3.5,-1.8))
lines(points_x_i,smooth_YX, col=rainbow(7)[1], lwd=2)
lines(points_x_i,smooth_YXC, col=rainbow(7)[2], lwd=2)
lines(points_x_i,smooth_YXU, col=rainbow(7)[3], lwd=2)
lines(points_x_i,smooth_YXCU, col=rainbow(7)[4], lwd=2)
#lines(points_x_i,smooth_YXCWZ, col=rainbow(7)[5], lwd=2)
lines(points_x_i,smooth_adj, col=rainbow(7)[6], lwd=2)
lines(points_x_i,smooth_adjC, col=rainbow(7)[7], lwd=2)
legend("topright",legend=c("YX","YXC","YXU","YXCU","adj","adjC"), col=rainbow(7)[-5], lwd=2)


###################################################################################

par(mfrow=c(1,1))
plot(data_analysis$X,data_analysis$Y, pch=19, cex=0.5,
     xlab="treatment: X", ylab="outcome: Y", ylim=c(-3.5,-1.8))
lines(points_x_i,smooth_YX, col=rainbow(7)[1], lwd=2)
lines(points_x_i,smooth_YXU, col=rainbow(7)[3], lwd=2)
lines(points_x_i,smooth_adj, col=rainbow(7)[6], lwd=2)
legend("topright",legend=c("YX","YXU","adj"), col=rainbow(7)[c(1,3,6)], lwd=2)

plot(data_analysis$X,data_analysis$Y, pch=19, cex=0.5,
     xlab="treatment: X", ylab="outcome: Y", ylim=c(-3.5,-1.8))
lines(points_x_i,smooth_YXC, col=rainbow(7)[2], lwd=2)
lines(points_x_i,smooth_YXCU, col=rainbow(7)[4], lwd=2)
lines(points_x_i,smooth_YXCWZ, col=rainbow(7)[5], lwd=2)
lines(points_x_i,smooth_adjC, col=rainbow(7)[7], lwd=2)
lines(points_x_i,smooth_adjCW, col=rainbow(7)[1], lwd=2)
legend("topright",legend=c("YXC","YXCU","YXCWZ","adjC","adjCW"), 
       col=rainbow(7)[c(2,4,5,7,1)], lwd=2)


plot(data_analysis$X,data_analysis$Y, pch=19, cex=0.5,
     xlab="treatment: X", ylab="outcome: Y", ylim=c(-3.5,-1.8))
lines(points_x_i,smooth_YXC, col=rainbow(7)[2], lwd=2)
lines(points_x_i,smooth_YXCU, col=rainbow(7)[4], lwd=2)
lines(points_x_i,smooth_adjC, col=rainbow(7)[7], lwd=2)
legend("topright",legend=c("YXC","YXCU","adjC"), col=rainbow(7)[c(2,4,7)], lwd=2)


###################################################################################

plot(data_analysis$X,data_analysis$Y, pch=19, cex=0.5,
     xlab="treatment: X", ylab="outcome: Y", 
     ylim=c(-3.5,-2.2), xlim=c(min_x+0.5,max_x-0.5))
for (i in seq(1,600,1)){
  lines(points_x_i,ksmooth(points_x_i,points_YX[i,],"normal", bandwidth = 0.5)$y,
        col="#F8710C7D", lwd=0.7, lty=1)
  lines(points_x_i,ksmooth(points_x_i,points_YXU[i,],"normal", bandwidth = 0.5)$y,
        col="#0BD97F7D", lwd=0.7, lty=1)
  lines(points_x_i,ksmooth(points_x_i,points_adj[i,],"normal", bandwidth = 0.5)$y,
        col="#0078EF7D", lwd=0.7, lty=1)
}
lines(points_x_i,smooth_YX, col=rainbow(7)[1], lwd=3)
lines(points_x_i,smooth_YXU, col=rainbow(7)[3], lwd=3)
lines(points_x_i,smooth_adj, col=rainbow(7)[6], lwd=3)
legend("topright",legend=c("YX","YXU","adj"), col=rainbow(7)[c(1,3,6)], lwd=2)


plot(data_analysis$X,data_analysis$Y, pch=19, cex=0.5,
      xlab="treatment: X", ylab="outcome: Y", ylim=c(-3.7,-1.5))
for (i in seq(1,1000,3)){
   lines(points_x_i,points_YX[i,], col=rainbow(7)[1], lwd=0.4, lty=1)
   lines(points_x_i,points_YXU[i,], col=rainbow(7)[3], lwd=0.4, lty=1)
   lines(points_x_i,points_adj[i,], col=rainbow(7)[6], lwd=0.4, lty=1)
}
legend("topright",legend=c("YX","YXU","adj"), col=rainbow(7)[c(1,3,6)], lwd=2)




plot(data_analysis$X,data_analysis$Y, pch=19, cex=0.5,
     xlab="treatment: X", ylab="outcome: Y", ylim=c(-3.5,-2))
for (i in seq(1,1000,3)){
  lines(points_x_i,ksmooth(points_x_i,points_YXU[i,],"normal", bandwidth = 0.5)$y,
        col=rainbow(7)[3], lwd=0.4, lty=1)
  #lines(points_x_i,points_YX[i,], col=rainbow(7)[3], lwd=0.4, lty=1)
}
for (i in seq(1,1000,3)){
  #lines(points_x_i,ksmooth(points_x_i,points_YX[i,],"normal", bandwidth = 0.15)$y,
  #      col=rainbow(7)[1], lwd=0.4, lty=1)
  lines(points_x_i,points_YX[i,], col=rainbow(7)[3], lwd=0.4, lty=1)
}


##########################################################################à


# --- intervals ---

par(mfrow=c(1,1))
plot(data_analysis$X,data_analysis$Y, pch=19, cex=0.4, col="gray",
     xlab="treatment: X", ylab="outcome: Y", axes=F,
     ylim=c(-3.5,-2.2), xlim=c(min_x,max_x))
for (i in seq(1,600,1)){
  #lines(points_x_i,ksmooth(points_x_i,points_YX[i,],"normal", bandwidth = 0.5)$y,
  #      col="#F8710C7D", lwd=0.7, lty=1)
  lines(points_x_i,ksmooth(points_x_i,points_YXU[i,],"normal", bandwidth = 0.5)$y,
        col="#0BD97F7D", lwd=0.7, lty=1)
  lines(points_x_i,ksmooth(points_x_i,points_adj[i,],"normal", bandwidth = 0.5)$y,
        col="#0078EF7D", lwd=0.7, lty=1)
}
#lines(points_x_i,smooth_YX, col=rainbow(7)[1], lwd=3)
lines(points_x_i,smooth_YXU, col=rainbow(7)[3], lwd=3)
lines(points_x_i,smooth_adj, col=rainbow(7)[6], lwd=3)
axis(side=2, at=c(-3.25,-3.00,-2.75,-2.50))
axis(side=1, at=c(7,9,11,13))
legend("topright",legend=c("YXU","adj"), col=rainbow(7)[c(3,6)], lwd=2,bty="n")


par(mfrow=c(1,1))
plot(data_analysis$X,data_analysis$Y, pch=19, cex=0.4, col="gray",
     xlab="treatment: X", ylab="outcome: Y", axes=F,
     ylim=c(-3.8,-2), xlim=c(min_x,max_x))
for (i in seq(1,600,2)){
  lines(points_x_i,ksmooth(points_x_i,points_YX[i,],"normal", bandwidth = 0.5)$y,
        col="#F8710C6D", lwd=1, lty=1)
  #lines(points_x_i,ksmooth(points_x_i,points_YXU[i,],"normal", bandwidth = 0.5)$y,
  #      col="#0BD97F7D", lwd=0.7, lty=1)
  lines(points_x_i,ksmooth(points_x_i,points_adj[i,],"normal", bandwidth = 0.5)$y,
        col="#0078EF6D", lwd=1, lty=1)
}
lines(points_x_i,smooth_YX, col=rainbow(7)[1], lwd=3)
#lines(points_x_i,smooth_YXU, col=rainbow(7)[3], lwd=3)
lines(points_x_i,smooth_adj, col=rainbow(7)[6], lwd=3)
axis(side=2, at=c(-3.25,-3.00,-2.75,-2.50))
axis(side=1, at=c(7,9,11,13))
legend("topright",legend=c("YX","adj"), col=rainbow(7)[c(1,6)], lwd=2,bty="n")



par(mfrow=c(1,1))
plot(data_analysis$X,data_analysis$Y, pch=19, cex=0.4, col="gray",
     xlab="treatment: X", ylab="outcome: Y", axes=F,
     ylim=c(-3.5,-2.2), xlim=c(min_x+1,max_x-0.5))
for (i in seq(1,600,1)){
  #lines(points_x_i,ksmooth(points_x_i,points_YXC[i,],"normal", bandwidth = 0.5)$y,
  #      col="#F8710C7D", lwd=0.7, lty=1)
  lines(points_x_i,ksmooth(points_x_i,points_YXCU[i,],"normal", bandwidth = 0.5)$y,
        col="#14DBD37D", lwd=0.7, lty=1)
  lines(points_x_i,ksmooth(points_x_i,points_adjC[i,],"normal", bandwidth = 0.5)$y,
        col="#EF58967D", lwd=0.7, lty=1)
}
#lines(points_x_i,smooth_YXC, col=rainbow(7)[2], lwd=2)
lines(points_x_i,smooth_YXCU, col=rainbow(7)[4], lwd=2)
lines(points_x_i,smooth_adjC, col=rainbow(7)[7], lwd=2)
axis(side=2, at=c(-3.25,-3.00,-2.75,-2.50))
axis(side=1, at=c(7,9,11,13))
legend("topright",legend=c("Y|XUC","adjC"), col=rainbow(7)[c(4,7)], lwd=2,bty="n")


par(mfrow=c(1,1))
plot(data_analysis$X,data_analysis$Y, pch=19, cex=0.4, col="gray",
     xlab="treatment: X", ylab="outcome: Y", axes=F,
     ylim=c(-3.5,-2.2), xlim=c(min_x+0.8,max_x-0.5))
for (i in seq(1,600,1)){
  lines(points_x_i,ksmooth(points_x_i,points_YXC[i,],"normal", bandwidth = 0.5)$y,
        col="#F8710C7D", lwd=0.7, lty=1)
  #lines(points_x_i,ksmooth(points_x_i,points_YXCU[i,],"normal", bandwidth = 0.5)$y,
  #      col="#14DBD37D", lwd=0.7, lty=1)
  lines(points_x_i,ksmooth(points_x_i,points_adjC[i,],"normal", bandwidth = 0.5)$y,
        col="#EF58967D", lwd=0.7, lty=1)
}
lines(points_x_i,smooth_YXC, col=rainbow(7)[2], lwd=2)
#lines(points_x_i,smooth_YXCU, col=rainbow(7)[4], lwd=2)
lines(points_x_i,smooth_adjC, col=rainbow(7)[7], lwd=2)
axis(side=2, at=c(-3.25,-3.00,-2.75,-2.50))
axis(side=1, at=c(7,9,11,13))
legend("topright",legend=c("Y|XC","adjC"), col=rainbow(7)[c(2,7)], lwd=2,bty="n")
