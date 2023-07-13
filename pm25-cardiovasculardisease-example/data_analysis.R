
source("BNP-For-Unmeasured-Confounding/src/models.R")
library(mvtnorm)
library(invgamma)
library(truncnorm)

load(
  "/n/dominici_nsaph_l3/Lab/projects/pm25-cardiovasculardisease-bnp/data_allariables_cleaned.RData"
)

data_clean$logPctblPvt2013 <- log(data_clean$PctblPvt2013)
data_clean$logPctblPvt2013[data_clean$logPctblPvt2013 == "-Inf"] = -4 
# add random noise to sum(data_clean$PctblPvt2013==0)

data_analysis = as.data.frame(
  cbind(
    Y = data_clean$lograte,
    X = data_clean$Exposure,
    U = data_clean$logMedHInc00_13,
    W = data_clean$PctOwnHs2013,
    Z = data_clean$logPctblPvt2013
  )
)

pdf("PctOwnHs2013-logMedHInc00_13-5quanitle.pdf")
cerf(data_analysis, X_split_method = split_5quantile, x_split_method = split_5quantile_x, n_group = 10, R = 2000, R_burnin = 1000, n=5362)
dev.off()
pdf("PctOwnHs2013-logMedHInc00_13-5quanitle2.pdf")
cerf(data_analysis, X_split_method = split_5quantile2, x_split_method = split_5quantile2_x, n_group = 10, R = 2000, R_burnin = 1000, n=5362)
dev.off()


data_analysis = as.data.frame(
  cbind(
    Y = data_clean$lograte,
    X = data_clean$Exposure,
    U = data_clean$logMedHInc00_13,
    W = data_clean$PctHighSchool,
    Z = data_clean$logPctblPvt2013
  )
)
pdf("PctHighSchool-logMedHInc00_13-5quanitle.pdf")
cerf(data_analysis, X_split_method = split_5quantile, x_split_method = split_5quantile_x, n_group = 10, R = 2000, R_burnin = 1000, n=5362)
dev.off()

pdf("PctHighSchool-logMedHInc00_13-5quanitle2.pdf")
cerf(data_analysis, X_split_method = split_5quantile2, x_split_method = split_5quantile2_x, n_group = 10, R = 2000, R_burnin = 1000, n=5362)
dev.off()


data_analysis = as.data.frame(
  cbind(
    Y = data_clean$lograte,
    X = data_clean$Exposure,
    U = data_clean$logMedHInc00_13,
    W = data_clean$PctOccupied,
    Z = data_clean$logPctblPvt2013
  )
)
pdf("PctOccupied-logPctblPvt2013-4quanitle.pdf")
cerf(data_analysis, X_split_method = split_4quantile, x_split_method = split_4quantile_x, n_group = 10, R = 2000, R_burnin = 1000, n=5362)
dev.off()

pdf("PctOccupied-logMedHInc00_13-5quanitle.pdf")
cerf(data_analysis, X_split_method = split_5quantile, x_split_method = split_5quantile_x, n_group = 10, R = 2000, R_burnin = 1000, n=5362)
dev.off()

pdf("PctOccupied-logMedHInc00_13-5quanitle2.pdf")
cerf(data_analysis, X_split_method = split_5quantile2, x_split_method = split_5quantile2_x, n_group = 10, R = 2000, R_burnin = 1000, n=5362)
dev.off()


data_analysis = as.data.frame(
  cbind(
    Y = data_clean$lograte,
    X = data_clean$Exposure,
    U = data_clean$logMedHInc00_13,
    W = data_clean$PctOccupied,
    Z = data_clean$logPctblPvt2013
  )
)

pdf("PctOccupied-logPctblPvt2013-4quanitle.pdf")
cerf(data_analysis, X_split_method = split_4quantile, x_split_method = split_4quantile_x, n_group = 10, R = 2000, R_burnin = 1000, n=5362)
dev.off()

pdf("PctOccupied-logPctblPvt2013-5quanitle.pdf")
cerf(data_analysis, X_split_method = split_5quantile, x_split_method = split_5quantile_x, n_group = 10, R = 2000, R_burnin = 1000, n=5362)
dev.off()


data_analysis = as.data.frame(
  cbind(
    Y = data_clean$lograte,
    X = data_clean$Exposure,
    U = data_clean$logMedHInc00_13,
    Z = data_clean$PctOccupied,
    W = data_clean$PctOwnHs2013
  )
)

pdf("PctOccupied-PctOwnHs2013-4quanitle.pdf")
cerf(data_analysis, X_split_method = split_4quantile, x_split_method = split_4quantile_x, n_group = 10, R = 2000, R_burnin = 1000, n=5362)
dev.off()

pdf("PctOccupied-PctOwnHs2013-5quanitle.pdf")
cerf(data_analysis, X_split_method = split_5quantile, x_split_method = split_5quantile_x, n_group = 10, R = 2000, R_burnin = 1000, n=5362)
dev.off()


pdf("PctOccupied-PctOwnHs2013-4quanitle2.pdf")
cerf(data_analysis, X_split_method = split_4quantile2, x_split_method = split_4quantile2_x, n_group = 10, R = 2000, R_burnin = 1000, n=5362)
dev.off()

pdf("PctOccupied-PctOwnHs2013-6quanitle.pdf")
cerf(data_analysis, X_split_method = split_6quantile, x_split_method = split_6quantile2_x, n_group = 10, R = 2000, R_burnin = 1000, n=5362)
dev.off()

data_analysis <- as.data.frame(cbind(
  Y = data_clean$lograte,
  X = data_clean$Exposure,
  U = data_clean$logMedHInc00_13,
  Z = data_clean$PctOccupied,
  W = data_clean$PctPoor
))


pdf("PctOccupied-PctPoor-4quanitle.pdf")
cerf(data_analysis, X_split_method = split_4quantile, x_split_method = split_4quantile_x, n_group = 10, R = 2000, R_burnin = 1000, n=5362)
dev.off()


pdf("PctOccupied-PctPoor-4quanitle2.pdf")
cerf(data_analysis, X_split_method = split_4quantile2, x_split_method = split_4quantile2_x, n_group = 10, R = 2000, R_burnin = 1000, n=5362)
dev.off()

pdf("PctOccupied-PctPoor-5quanitle.pdf")
cerf(data_analysis, X_split_method = split_5quantile, x_split_method = split_5quantile_x, n_group = 10, R = 1000, R_burnin = 500, n=5362)
dev.off()
###################################################################################

# --- estimation: ---



cerf <- function(data_analysis,
                 X_split_method,
                 x_split_method,
                 n_group = 10,
                 # max number of groups
                 R = 1000,
                 # iterations
                 R_burnin = 500,
                 # burn-in
                 n = 5362) {
  # model Y|X
  X_tilde_YX = cbind(rep(1, n), data_analysis$X)
  
  start.time <- Sys.time()
  model_YX = DDP_NOadj(
    s_seed = 2,
    X_tilde_YX,
    X_split_method = X_split_method,
    R = R,
    R_burnin = R_burnin,
    n_group = n_group,
    n = n,
    data_analysis = data_analysis
  )
  Sys.time() - start.time
  
  # model Y|X,U
  X_tilde_YXU = cbind(rep(1, n), data_analysis$X, data_analysis$U)
  X_tilde_YXU = data.matrix(X_tilde_YXU)
  
  start.time <- Sys.time()
  model_YXU = DDP_NOadj(
    s_seed = 1,
    X_tilde_YXU,
    X_split_method = X_split_method,
    R = R,
    R_burnin = R_burnin,
    n_group = n_group,
    n = n,
    data_analysis = data_analysis
  )
  Sys.time() - start.time
  
  
  # model adjusted: Y|X,Z and W|X,Z
  X_tilde_adj = cbind(rep(1, n), data_analysis$X, data_analysis$Z)
  X_tilde_adj = data.matrix(X_tilde_adj)
  
  start.time <- Sys.time()
  model_adj = DDP_ADJ(
    s_seed = 1,
    X_tilde_adj,
    X_tilde_adj,
    X_split_method = X_split_method,
    R = R,
    R_burnin = R_burnin,
    n_group = n_group,
    n = n,
    data_analysis = data_analysis
  )
  Sys.time() - start.time
  
  
  min_x = min(data_analysis$X)
  max_x = max(data_analysis$X)
  
  # causal effect for a grid of points
  points_x_i = seq(min_x, max_x, length.out = 100)
  
  points_YX = sapply(points_x_i, function(x)
    curve_NOadj(
      x,
      post_chain = model_YX,
      X_tilde = X_tilde_YX,
       x_split_method = x_split_method,
     data_analysis = data_analysis,
      n_group = n_group 
      ))
  points_YXU = sapply(points_x_i, function(x)
    curve_NOadj(
      x,
      post_chain = model_YXU,
      X_tilde = X_tilde_YXU,
      x_split_method = x_split_method,
     data_analysis = data_analysis,
      n_group = n_group 
    ))
  
  points_adj = sapply(points_x_i, function(x)
    curve_ADJ(
      x,
      post_chain = model_adj,
      X_tilde = X_tilde_adj,
       x_split_method = x_split_method,
       data_analysis = data_analysis,
       n_group = n_group 
    ))
  
  smooth_YX = ksmooth(
    points_x_i,
    apply(points_YX, 2, quantile, prob = 0.5, na.rm = TRUE),
    "normal",
    bandwidth = 0.5
  )$y
  
  smooth_YXU = ksmooth(
    points_x_i,
    apply(points_YXU, 2, quantile, prob = 0.5, na.rm =
            TRUE),
    "normal",
    bandwidth = 0.5
  )$y
  
  smooth_adj = ksmooth(
    points_x_i,
    apply(points_adj, 2, quantile, prob = 0.5, na.rm =
            TRUE),
    "normal",
    bandwidth = 0.5
  )$y
  
  par(mfrow = c(1, 1))
  plot(
    data_analysis$X,
    data_analysis$Y,
    pch = 19,
    cex = 0.3,
    col = 'dark grey',
    xlab = "treatment: X",
    ylab = "outcome: Y",
    ylim = c(-3.5,-1.8)
  )
  lines(points_x_i, smooth_YX, col = rainbow(7)[1], lwd = 2)
  lines(points_x_i, smooth_YXU, col = rainbow(7)[3], lwd = 2)
  lines(points_x_i, smooth_adj, col = rainbow(7)[6], lwd = 2)
  legend(
    "topright",
    legend = c("YX", "YXU", "adj"),
    col = rainbow(7)[c(1, 3, 6)],
    lwd = 2
  )
}