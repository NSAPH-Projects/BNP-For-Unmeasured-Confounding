###################################################################################
# ---   Compare Estimated Causal Exposure-Response Function (CERF)           ---
# ---   Compare 1) adjsut for confounder U directly                          ---
# ---           2) not adjsut for confounder U at all                        ---
# ---           3) adjust for confounder U with negative controls labeled Z and 
#                  W when U is masked                                          ---
###################################################################################

source("src/BNP_NC.R")
source("src/BNP_NoAdj.R")
source("src/helper.R")

cerf <- function(data, 
                 # a dataframe including columns of exposure X, outcome Y, 
                 #  confounder U, negative control exposure Z, and negative 
                 #  control outcome W
                 X_split_method = "split_4quantile",
                 # how we divide the training data of X when fitting the model
                 x_split_method = "split_4quantile_x",
                 # how we divide the testing data of x when applying the model
                 probs = NULL, 
                 # a vector of values between [0,1], needed when split_quantile is
                 # called as x_split method
                 pts = NULL, 
                 # a vector of values, needed when split_at_fixed_pt is called as 
                 # x_split method
                 n_group = 10,
                 # max number of groups
                 R = 1000,
                 # number of iterations
                 R_burnin = 500,
                 # burn-in
                 n, 
                 # sample size,
                 s_seed = 3
                 # seed number
                 ) {
  # model Y|X
  X_tilde_YX = cbind(rep(1, n), data$X)
  
  start.time <- Sys.time()
  model_YX = DDP_NOadj(
    s_seed = s_seed,
    X_tilde_YX,
    X_split_method = X_split_method,
    probs = probs, 
    pts = pts, 
    R = R,
    R_burnin = R_burnin,
    n_group = n_group,
    n = n,
    data = data
  )
  Sys.time() - start.time
  
  # model Y|X,U
  X_tilde_YXU = cbind(rep(1, n), data$X, data$U)
  X_tilde_YXU = data.matrix(X_tilde_YXU)
  
  start.time <- Sys.time()
  model_YXU = DDP_NOadj(
    s_seed = s_seed, 
    X_tilde_YXU,
    X_split_method = X_split_method,
    probs = probs, 
    pts = pts, 
    R = R,
    R_burnin = R_burnin,
    n_group = n_group,
    n = n,
    data = data
  )
  Sys.time() - start.time
  
  
  # model Y|X, Z, W
  X_tilde_adj = cbind(rep(1, n), data$X, data$Z)
  X_tilde_adj = data.matrix(X_tilde_adj)
  
  start.time <- Sys.time()
  model_adj = DDP_ADJ(
    s_seed = s_seed,
    X_tilde_adj,
    X_tilde_adj,
    X_split_method = X_split_method,
    probs = probs, 
    pts = pts, 
    R = R,
    R_burnin = R_burnin,
    n_group = n_group,
    n = n,
    data = data
  )
  Sys.time() - start.time
  
  
  # obtain the estimated causal effect for a grid of points of x using 
  # the three fitted models
  
  min_x = min(data$X)
  max_x = max(data$X)
  points_x_i = seq(min_x, max_x, length.out = 100)
  
  # based on model Y|X
  points_YX = sapply(points_x_i, function(x)
    curve_NOadj(
      x,
      post_chain = model_YX,
      X_tilde = X_tilde_YX,
      x_split_method = x_split_method,
      probs = probs, 
      pts = pts, 
      data = data,
      n_group = n_group 
    ))
  
  smooth_YX = ksmooth(
    points_x_i,
    apply(points_YX, 2, quantile, prob = 0.5, na.rm = TRUE),
    "normal",
    bandwidth = 0.5
  )$y
  
  # based on model Y|X, U
  points_YXU = sapply(points_x_i, function(x)
    curve_NOadj(
      x,
      post_chain = model_YXU,
      X_tilde = X_tilde_YXU,
      x_split_method = x_split_method,
      probs = probs, 
      pts = pts, 
      data = data,
      n_group = n_group 
    ))

  smooth_YXU = ksmooth(
    points_x_i,
    apply(points_YXU, 2, quantile, prob = 0.5, na.rm =
            TRUE),
    "normal",
    bandwidth = 0.5
  )$y
  
  # based on model Y|X, Z,W
  points_adj = sapply(points_x_i, function(x)
    curve_ADJ(
      x,
      post_chain = model_adj,
      X_tilde = X_tilde_adj,
      x_split_method = x_split_method,
      probs = probs, 
      pts = pts, 
      data = data,
      n_group = n_group 
    ))
  
  smooth_adj = ksmooth(
    points_x_i,
    apply(points_adj, 2, quantile, prob = 0.5, na.rm =
            TRUE),
    "normal",
    bandwidth = 0.5
  )$y
  return(list(points_YX=points_YX,
              points_YXU=points_YXU,
              points_adj=points_adj,
              points_x_i=points_x_i,
              smooth_YX=smooth_YX,
              smooth_YXU=smooth_YXU,
              smooth_adj=smooth_adj))
}