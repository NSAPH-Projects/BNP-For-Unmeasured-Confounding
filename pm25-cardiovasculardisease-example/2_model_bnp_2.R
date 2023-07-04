###################################################################################
# ---          DATA ANALYSIS FOR THE PAPER         ---
# ---                   4 settings                 ---
###################################################################################

# library
library(mvtnorm)
library(invgamma)
library(BNPmix)
library(truncnorm)
library(ggplot2)
library(GGally)
library(corrgram)

###################################################################################

# uploasd dataset
#load("/n/dominici_nsaph_l3/Lab/projects/pm25-cardiovasculardisease-bnp/data_cleaned.RData")
load(
  "/n/dominici_nsaph_l3/Lab/projects/pm25-cardiovasculardisease-bnp/data_allariables_cleaned.RData"
)

###################################################################################

# cut the extremes for the treatment
summary(data_clean$Exposure)
plot(data_clean$Exposure, data_clean$lograte)
abline(v = 5, col = "red")
abline(v = 14, col = "red")

data_clean <-
  data_clean[data_clean$Exposure > 5 & data_clean$Exposure < 14, ]

# checking correlation among other variables
corrgram(
  data_clean[, c(4:28)],
  order = TRUE,
  lower.panel = panel.shade,
  upper.panel = panel.conf
)
covariates <- c("PctFemale",
                "smokerate2013",
                "AvgCommute",
                "BMI2013",
                "avgASOStemp",
                "mean_age")
corrgram(
  data_clean[, covariates],
  order = TRUE,
  lower.panel = panel.shade,
  upper.panel = panel.conf
)

# preparing dataset for the Gibbs

data_analysis = as.data.frame(
  cbind(
    Y = data_clean$lograte,
    X = data_clean$Exposure,
    U = data_clean$logMedHInc00_13,
    W = data_clean$PctHighSchool,
    Z = exp(data_clean$PctWhite),
    data_clean[, covariates]
  )
)


# correlation
#pdf("corr.pdf",width=8, height=8)
#ggpairs(data_analysis)
#dev.off()

###################################################################################
# ---    GIBBS for Dependent Dirichlet Mixture Model      ----
###################################################################################

# basic setting for Gibbs

R = 2000            # iteartions
R_burnin = 1000      # burn-in
n_group = 10        # max number of groups

n = length(data_analysis$Y)

###################################################################################

# General functions for clustering

split_4quantile<-function(x_matrix){
  seq_X=c(min(x_matrix),quantile(x_matrix,c(0.25,0.5,0.75)), max(x_matrix))
  X_weights=cbind(x_matrix*I(x_matrix<seq_X[2]),
                  x_matrix*I(x_matrix<seq_X[3] & x_matrix>seq_X[2]),
                  x_matrix*I(x_matrix<seq_X[4] & x_matrix>seq_X[3]),
                  x_matrix*I(x_matrix>seq_X[4]))
  X_weights=(X_weights!=0)*(X_weights)+
    (X_weights==0)*matrix(rnorm(prod(dim(X_weights)),0,0.003),
                          ncol=dim(X_weights)[2],nrow=dim(X_weights)[1])
  
  return(X_weights)
}

###################################################################################

# Gibbs samplers

DDP_NOadj <- function(s_seed = 1, X_tilde) {
  set.seed(s_seed)
  dim_REG = dim(X_tilde)[2]
  
  #prior
  alpha_mu = rep(0, dim_REG)
  alpha_sigma = 1
  gamma_y_prior = c(1, 0.5)
  eta_prior = c(0, 5)
  
  #initialization
  alpha = matrix(rep(rmvnorm(
    1, alpha_mu * rep(1, dim_REG),
    alpha_sigma * diag(dim_REG)
  ), n_group),
  ncol = n_group, nrow = dim_REG)
  gamma_y = rep(rinvgamma(1, gamma_y_prior[1], gamma_y_prior[2]), n_group)
  eta = matrix(
    rnorm((n_group - 1) * 2, eta_prior[1], eta_prior[2]),
    ncol = n_group - 1,
    nrow = dim(X_weights)[2]
  )
  eta_X = pnorm(cbind(X_weights %*% (eta[, 1]), X_weights %*% (eta[, 2])))
  for (g_k in 3:(n_group - 1)) {
    eta_X = cbind(eta_X, pnorm(X_weights %*% (eta[, g_k])))
  }
  eta_X = cbind(eta_X, 1)
  tau_prov = cbind(
    eta_X[, 1] * dnorm(data_analysis$Y, X_tilde %*% (alpha[, 1]), sqrt(gamma_y[1])),
    eta_X[, 2] * (1 - eta_X[, 1]) * dnorm(data_analysis$Y, X_tilde %*%
                                            (alpha[, 2]), sqrt(gamma_y[2]))
  )
  for (g_k in 3:n_group) {
    tau_prov = cbind(
      tau_prov,
      eta_X[, g_k] * apply(1 - eta_X[, 1:(g_k - 1)], 1, prod) *
        dnorm(data_analysis$Y, X_tilde %*% (alpha[, g_k]), sqrt(gamma_y[g_k]))
    )
  }
  tau = tau_prov / apply(tau_prov, 1, sum)
  
  K = sample(1:n_group, n, replace = TRUE)
  table_Ki = rep(0, n_group)
  table_Ki[sort(unique(K))] = table(K)
  
  Q = matrix(NA, nrow = n, ncol = n_group - 1)
  
  #chains
  post_alpha = matrix(NA, nrow = n_group * dim_REG, ncol = R)
  post_eta = matrix(NA, nrow = (n_group - 1) * (dim(X_weights)[2]), ncol =
                      R)
  post_beta_x = matrix(NA, nrow = n_group, ncol = R)
  post_intercept = matrix(NA, nrow = n_group, ncol = R)
  
  for (r in 1:R) {
    # ----- model for Y given variables : ----
    
    #mu_Q
    mu_Q = cbind((tau[, 1]), (tau[, 2] / (1 - tau[, 1])))
    if (n_group > 3) {
      mu_Q = cbind(mu_Q, sapply(3:(n_group - 1), function(l)
        (tau[, l] / (
          1 - apply(tau[, 1:(l - 1)], 1, sum)
        ))))
    }
    # for resolve problems in the bounds
    mu_Q[which(is.nan(mu_Q))] = 1
    mu_Q[which(mu_Q > 1)] = 1
    mu_Q = mu_Q - 9.9e-15 * (mu_Q > (1 - 1e-16))
    mu_Q = mu_Q + 1e-16 * (mu_Q < 1e-16)
    
    #Q
    for (i in 1:n) {
      for (l in 1:(min(K[i], n_group - 1))) {
        if (l < K[i]) {
          Q[i, l] = rtruncnorm(1, b = 0, mean = qnorm(mu_Q[i, l]))
        } else{
          Q[i, l] = rtruncnorm(1, a = 0, mean = qnorm(mu_Q[i, l]))
        }
      }
    }
    
    #eta
    check_q = sapply(1:(n_group - 1), function(l)
      length(Q[K >= l, l]))
    v_eta = lapply(which(check_q > 2), function(l)
      1 / eta_prior[2] + t(X_weights[K >= l, ]) %*% X_weights[K >= l, ])
    m_eta = sapply(which(check_q > 2), function(l)
      eta_prior[1] / eta_prior[2] +
        t(X_weights[K >= l, ]) %*% Q[K >= l, l] / gamma_y[l])
    
    if (any(check_q < 3)) {
      for (l in which(check_q == 1 | check_q == 2)) {
        v_eta[[l]] = diag(dim(X_weights)[2]) / eta_prior[2]
        m_eta = cbind(m_eta, rep(eta_prior[1] / eta_prior[2], 2))
      }
    }
    eta[, which(table_Ki[-n_group] > 0)] = sapply(which(table_Ki[-n_group] >
                                                          0),
                                                  function(l)
                                                    rmvnorm(1, solve(v_eta[[l]]) %*% m_eta[, l]))
    eta[, which(table_Ki[-n_group] == 0)] = rnorm(length(which(table_Ki[-n_group] ==
                                                                 0)) * 2, eta_prior[1], sqrt(eta_prior[2]))
    
    
    #tau ( recursive weigths in the mixture )
    eta_X = pnorm(cbind(X_weights %*% (eta[, 1]), X_weights %*% (eta[, 2])))
    for (g_k in 3:(n_group - 1)) {
      eta_X = cbind(eta_X, pnorm(X_weights %*% (eta[, g_k])))
    }
    eta_X = cbind(eta_X, 1)
    
    tau_prov = cbind(
      eta_X[, 1] * dnorm(data_analysis$Y, X_tilde %*% (alpha[, 1]), sqrt(gamma_y[1])),
      eta_X[, 2] * (1 - eta_X[, 1]) * dnorm(data_analysis$Y, X_tilde %*%
                                              (alpha[, 2]), sqrt(gamma_y[2]))
    )
    for (g_k in 3:n_group) {
      tau_prov = cbind(
        tau_prov,
        eta_X[, g_k] * apply(1 - eta_X[, 1:(g_k - 1)], 1, prod) *
          dnorm(data_analysis$Y, X_tilde %*% (alpha[, g_k]), sqrt(gamma_y[g_k]))
      )
    }
    
    tau = tau_prov / apply(tau_prov, 1, sum)
    for (i in which(is.na(tau[, 1]))) {
      tau[i, ] = 1 / n_group
    }
    
    #K
    K = sapply(1:n, function(i)
      (1:n_group) %*% rmultinom(1, 1, tau[i, ]))
    table_Ki = rep(0, n_group)
    table_Ki[sort(unique(K))] = table(K)
    
    #alpha + gamma_y ( parameters of the mean + variance in each cluster of the distribution of Y|-- )
    for (g_k in 1:n_group) {
      if (length(which(K == g_k)) == 1) {
        V = diag(dim_REG) / alpha_sigma + t(t(X_tilde[K == g_k, ])) %*% t(X_tilde[K ==
                                                                                    g_k, ]) / gamma_y[g_k]
        M = alpha_mu / alpha_sigma + t(t(X_tilde[K == g_k, ])) %*% t(data_analysis$Y[K ==
                                                                                       g_k]) / gamma_y[g_k]
      } else{
        V = diag(dim_REG) / alpha_sigma + t(X_tilde[K == g_k, ]) %*% X_tilde[K ==
                                                                               g_k, ] / gamma_y[g_k]
        M = alpha_mu / alpha_sigma + t(X_tilde[K == g_k, ]) %*% data_analysis$Y[K ==
                                                                                  g_k] / gamma_y[g_k]
      }
      alpha[, g_k] = t(rmvnorm(1, solve(V) %*% M, solve(V)))
      
      G = sum((data_analysis$Y[K == g_k] - X_tilde[K == g_k, ] %*% (alpha[, g_k])) ^
                2)
      gamma_y[g_k] = rinvgamma(1, gamma_y_prior[1] + (sum(K == g_k)) / 2, gamma_y_prior[2] +
                                 G / 2)
    }
    
    
    #beta_x + intercept
    beta_x = alpha[2, ]
    if (dim_REG < 3) {
      intercept = sapply(1:n_group, function(g)
        alpha[1, g])
    } else{
      intercept = sapply(1:n_group, function(g)
        alpha[1, g] +
          apply(t(t(X_tilde[, -c(1, 2)])), 2, mean) %*% (alpha[-c(1, 2), g]))
    }
    
    
    # ---------- save the chains  ------
    post_alpha[, r] = c(alpha)
    post_eta[, r] = c(eta)
    post_beta_x[, r] = beta_x
    post_intercept[, r] = intercept
    
    #check
    if (r %% 500 == 0)
      print(r)
  }
  
  print(paste0("sample ", s_seed))
  
  return(
    list(
      post_alpha = post_alpha[, (R_burnin + 1):R],
      post_beta_x = post_beta_x[, (R_burnin + 1):R],
      post_intercept = post_intercept[, (R_burnin + 1):R],
      post_eta = post_eta[, (R_burnin + 1):R]
    )
  )
}

DDP_ADJ <- function(s_seed = 1, X_tilde, X_W_reg) {
  set.seed(s_seed)
  dim_REG = dim(X_tilde)[2]
  dim_W = dim(X_W_reg)[2]
  
  #prior
  alpha_mu = rep(0, dim_REG)
  alpha_sigma = 1
  delta_mu = rep(0, dim_W)
  delta_sigma = 1
  gamma_y_prior = c(1, 0.5)
  gamma_w_prior = c(2, 2)
  eta_prior = c(0, 5)
  
  #initialization
  alpha = matrix(rep(rmvnorm(
    1, alpha_mu * rep(1, dim_REG),
    alpha_sigma * diag(dim_REG)
  ), n_group),
  ncol = n_group, nrow = dim_REG)
  delta = rnorm(dim_W, delta_mu, delta_sigma)
  gamma_y = rep(rinvgamma(1, gamma_y_prior[1], gamma_y_prior[2]), n_group)
  gamma_w = rinvgamma(1, gamma_w_prior[1], gamma_w_prior[2])
  eta = matrix(
    rnorm((n_group - 1) * 2, eta_prior[1], eta_prior[2]),
    ncol = n_group - 1,
    nrow = dim(X_weights)[2]
  )
  eta_X = pnorm(cbind(X_weights %*% (eta[, 1]), X_weights %*% (eta[, 2])))
  for (g_k in 3:(n_group - 1)) {
    eta_X = cbind(eta_X, pnorm(X_weights %*% (eta[, g_k])))
  }
  eta_X = cbind(eta_X, 1)
  tau_prov = cbind(
    eta_X[, 1] * dnorm(data_analysis$Y, X_tilde %*% (alpha[, 1]), sqrt(gamma_y[1])),
    eta_X[, 2] * (1 - eta_X[, 1]) * dnorm(data_analysis$Y, X_tilde %*%
                                            (alpha[, 2]), sqrt(gamma_y[2]))
  )
  for (g_k in 3:n_group) {
    tau_prov = cbind(
      tau_prov,
      eta_X[, g_k] * apply(1 - eta_X[, 1:(g_k - 1)], 1, prod) *
        dnorm(data_analysis$Y, X_tilde %*% (alpha[, g_k]), sqrt(gamma_y[g_k]))
    )
  }
  tau = tau_prov / apply(tau_prov, 1, sum)
  
  K = sample(1:n_group, n, replace = TRUE)
  table_Ki = rep(0, n_group)
  table_Ki[sort(unique(K))] = table(K)
  
  Q = matrix(NA, nrow = n, ncol = n_group - 1)
  beta_x = rep(0, n_group)
  
  #chains
  post_alpha = matrix(NA, nrow = n_group * dim_REG, ncol = R)
  post_delta = matrix(NA, nrow = dim_W, ncol = R)
  post_beta_x = matrix(NA, nrow = n_group, ncol = R)
  post_intercept = matrix(NA, nrow = n_group, ncol = R)
  post_eta = matrix(NA, nrow = (n_group - 1) * (dim(X_weights)[2]), ncol =
                      R)
  
  for (r in 1:R) {
    # ----- model for W : ----
    
    #delta  ( parameters of the mean in the distribution of W|X,Z )
    V_delta = diag(dim_W) / delta_sigma + t(X_W_reg) %*% X_W_reg / gamma_w
    m_delta = delta_mu / delta_sigma + t(X_W_reg) %*% data_analysis$W /
      gamma_w
    delta = rmvnorm(1, solve(V_delta) %*% m_delta, solve(V_delta))
    
    #gamma_w    ( parameter of the variance in the distribution of W|X,Z )
    g_w = sum((data_analysis$W - X_W_reg %*% t(delta)) ^ 2)
    gamma_w = rinvgamma(1, gamma_w_prior[1] + n / 2, gamma_w_prior[2] +
                          g_w / 2)
    
    #for compute eta (regression parameters in the weights)
    # we need to compute mu_Q and Q
    
    # ----- model for Y with proxies : ----
    
    #mu_Q
    mu_Q = cbind((tau[, 1]), (tau[, 2] / (1 - tau[, 1])))
    if (n_group > 3) {
      mu_Q = cbind(mu_Q, sapply(3:(n_group - 1), function(l)
        (tau[, l] / (
          1 - apply(tau[, 1:(l - 1)], 1, sum)
        ))))
    }
    # for resolve problems in the bounds
    mu_Q[which(is.nan(mu_Q))] = 1
    mu_Q[which(mu_Q > 1)] = 1
    mu_Q = mu_Q - 9.9e-15 * (mu_Q > (1 - 1e-16))
    mu_Q = mu_Q + 1e-16 * (mu_Q < 1e-16)
    
    #Q
    for (i in 1:n) {
      for (l in 1:(min(K[i], n_group - 1))) {
        if (l < K[i]) {
          Q[i, l] = rtruncnorm(1, b = 0, mean = qnorm(mu_Q[i, l]))
        } else{
          Q[i, l] = rtruncnorm(1, a = 0, mean = qnorm(mu_Q[i, l]))
        }
      }
    }
    
    #eta
    check_q = sapply(1:(n_group - 1), function(l)
      length(Q[K >= l, l]))
    v_eta = lapply(which(check_q > 2), function(l)
      1 / eta_prior[2] + t(X_weights[K >= l, ]) %*% X_weights[K >= l, ])
    m_eta = sapply(which(check_q > 2), function(l)
      eta_prior[1] / eta_prior[2] +
        t(X_weights[K >= l, ]) %*% Q[K >= l, l] / gamma_y[l])
    
    if (any(check_q < 3)) {
      for (l in which(check_q == 1 | check_q == 2)) {
        v_eta[[l]] = diag(dim(X_weights)[2]) / eta_prior[2]
        m_eta = cbind(m_eta, rep(eta_prior[1] / eta_prior[2], 2))
      }
    }
    eta[, which(table_Ki[-n_group] > 0)] = sapply(which(table_Ki[-n_group] >
                                                          0),
                                                  function(l)
                                                    rmvnorm(1, solve(v_eta[[l]]) %*% m_eta[, l]))
    eta[, which(table_Ki[-n_group] == 0)] = rnorm(length(which(table_Ki[-n_group] ==
                                                                 0)) * 2, eta_prior[1], sqrt(eta_prior[2]))
    
    
    #tau ( recursive weigths in the mixture )
    eta_X = pnorm(cbind(X_weights %*% (eta[, 1]), X_weights %*% (eta[, 2])))
    for (g_k in 3:(n_group - 1)) {
      eta_X = cbind(eta_X, pnorm(X_weights %*% (eta[, g_k])))
    }
    eta_X = cbind(eta_X, 1)
    
    tau_prov = cbind(
      eta_X[, 1] * dnorm(data_analysis$Y, X_tilde %*% (alpha[, 1]), sqrt(gamma_y[1])),
      eta_X[, 2] * (1 - eta_X[, 1]) * dnorm(data_analysis$Y, X_tilde %*%
                                              (alpha[, 2]), sqrt(gamma_y[2]))
    )
    for (g_k in 3:n_group) {
      tau_prov = cbind(
        tau_prov,
        eta_X[, g_k] * apply(1 - eta_X[, 1:(g_k - 1)], 1, prod) *
          dnorm(data_analysis$Y, X_tilde %*% (alpha[, g_k]), sqrt(gamma_y[g_k]))
      )
    }
    
    tau = tau_prov / apply(tau_prov, 1, sum)
    for (i in which(is.na(tau[, 1]))) {
      tau[i, ] = 1 / n_group
    }
    
    #K
    K = sapply(1:n, function(i)
      (1:n_group) %*% rmultinom(1, 1, tau[i, ]))
    table_Ki = rep(0, n_group)
    table_Ki[sort(unique(K))] = table(K)
    
    #alpha + gamma_y ( parameters of the mean + variance in each cluster of the distribution of Y|X,Z )
    for (g_k in 1:n_group) {
      if (length(which(K == g_k)) == 1) {
        V = diag(dim_REG) / alpha_sigma + t(t(X_tilde[K == g_k, ])) %*% t(X_tilde[K ==
                                                                                    g_k, ]) / gamma_y[g_k]
        M = alpha_mu / alpha_sigma + t(t(X_tilde[K == g_k, ])) %*% t(data_analysis$Y[K ==
                                                                                       g_k]) / gamma_y[g_k]
      } else{
        V = diag(dim_REG) / alpha_sigma + t(X_tilde[K == g_k, ]) %*% X_tilde[K ==
                                                                               g_k, ] / gamma_y[g_k]
        M = alpha_mu / alpha_sigma + t(X_tilde[K == g_k, ]) %*% data_analysis$Y[K ==
                                                                                  g_k] / gamma_y[g_k]
      }
      alpha[, g_k] = t(rmvnorm(1, solve(V) %*% M, solve(V)))
      
      G = sum((data_analysis$Y[K == g_k] - X_tilde[K == g_k, ] %*% (alpha[, g_k])) ^
                2)
      gamma_y[g_k] = rinvgamma(1, gamma_y_prior[1] + (sum(K == g_k)) / 2, gamma_y_prior[2] +
                                 G / 2)
    }
    
    
    #beta_x + intercept
    beta_x = alpha[2, ] - alpha[3, ] * delta[1, 2] / delta[1, 3]
    if (dim_REG < 4) {
      intercept = sapply(1:n_group, function(g)
        mean(
          alpha[1, g] + alpha[3, g] * mean(data_analysis$Z) +
            alpha[3, g] * delta[1, 2] / delta[1, 3] * mean(data_analysis$X)
        ))
    } else{
      intercept = sapply(1:n_group, function(g)
        mean(
          alpha[1, g] + alpha[3, g] * mean(data_analysis$Z) +
            alpha[3, g] * delta[1, 2] / delta[1, 3] * mean(data_analysis$X)
        ) +
          sum(alpha[4:dim_REG, g] * apply(X_tilde[, 4:dim_REG], 2, mean)))
    }
    
    # ---------- save the chains  ------
    post_alpha[, r] = c(alpha)
    post_delta[, r] = delta[, 1:dim_W]
    post_beta_x[, r] = beta_x
    post_intercept[, r] = intercept
    post_eta[, r] = c(eta)
    
    #check
    if (r %% 500 == 0)
      print(r)
  }
  
  print(paste0("sample ", s_seed))
  
  return(
    list(
      post_alpha = post_alpha[, (R_burnin + 1):R],
      post_beta_x = post_beta_x[, (R_burnin + 1):R],
      post_intercept = post_intercept[, (R_burnin + 1):R],
      post_eta = post_eta[, (R_burnin + 1):R],
      post_delta = post_delta[, (R_burnin + 1):R]
    )
  )
}

###################################################################################

# --- estimation: ---

# common part
X_weights = split_4quantile(data_analysis$X)


# model Y|X
X_tilde_YX = cbind(rep(1, n), data_analysis$X)

start.time <- Sys.time()
model_YX = DDP_NOadj(s_seed = 1, X_tilde_YX)
Sys.time() - start.time

# model Y|X,U
X_tilde_YXU = cbind(rep(1, n), data_analysis$X, data_analysis$U)
X_tilde_YXU = data.matrix(X_tilde_YXU)

start.time <- Sys.time()
model_YXU = DDP_NOadj(s_seed = 1, X_tilde_YXU)
Sys.time() - start.time


# model adjusted: Y|X,Z and W|X,Z
X_tilde_adj = cbind(rep(1, n), data_analysis$X, data_analysis$Z)
X_tilde_adj = data.matrix(X_tilde_adj)

start.time <- Sys.time()
model_adj = DDP_ADJ(s_seed = 1, X_tilde_adj, X_tilde_adj)
Sys.time() - start.time


# model Y|X,cov
X_tilde_YXC = cbind(rep(1, n), data_analysis$X, data_analysis[, covariates])
X_tilde_YXC = data.matrix(X_tilde_YXC)

start.time <- Sys.time()
model_YXC = DDP_NOadj(s_seed = 1, X_tilde_YXC)
Sys.time() - start.time

# model Y|X,cov,U
X_tilde_YXCU = cbind(rep(1, n), data_analysis$X, data_analysis[, covariates],
                     data_analysis$U)
X_tilde_YXCU = data.matrix(X_tilde_YXCU)

start.time <- Sys.time()
model_YXCU = DDP_NOadj(s_seed = 1, X_tilde_YXCU)
Sys.time() - start.time


# model Y|X,cov,W,Z
X_tilde_YXCWZ = cbind(rep(1, n),
                      data_analysis$X,
                      data_analysis[, covariates],
                      data_analysis$W,
                      data_analysis$Z)
X_tilde_YXCWZ = data.matrix(X_tilde_YXCWZ)

start.time <- Sys.time()
model_YXCWZ = DDP_NOadj(s_seed = 1, X_tilde_YXCWZ)
Sys.time() - start.time


# model adjusted: Y|X,Z,cov and W|X,Z
X_tilde_adjC = cbind(rep(1, n), data_analysis$X, data_analysis$Z,
                     data_analysis[, covariates])
X_tilde_adjC = data.matrix(X_tilde_adjC)
X_W = data.matrix(cbind(rep(1, n), data_analysis$X, data_analysis$Z))

start.time <- Sys.time()
model_adjC = DDP_ADJ(s_seed = 2, X_tilde_adjC, X_W)
Sys.time() - start.time

# model adjusted: Y|X,Z,cov and W|X,Z,cov
start.time <- Sys.time()
model_adjCW = DDP_ADJ(s_seed = 2, X_tilde_adjC, X_tilde_adjC)
Sys.time() - start.time


#save.image("results_models_data_2.RData")

###################################################################################
# ---    compute MEDIANS for the posterior distributions      ----
###################################################################################


split_4quantile_x <- function(x_matrix, x) {
  seq_X = c(min(x_matrix), quantile(x_matrix, c(0.25, 0.5, 0.75)), max(x_matrix))
  X_weights = cbind(x * I(x < seq_X[2]),
                    x * I(x < seq_X[3] & x > seq_X[2]),
                    x * I(x < seq_X[4] & x > seq_X[3]),
                    x * I(x > seq_X[4]))
  X_weights = (X_weights != 0) * (X_weights) +
    (X_weights == 0) * matrix(rnorm(prod(dim(X_weights)), 0, 0.003),
                              ncol = dim(X_weights)[2],
                              nrow = dim(X_weights)[1])
  
  return(X_weights)
}

#################################################################################

curve_NOadj <- function(x, post_chain, X_tilde) {
  X_weights = split_4quantile_x(data_analysis$X, x)
  dim_Xw = dim(X_weights)[2]
  dim_REG = dim(X_tilde)[2]
  
  c_eta_X = pnorm(cbind(t(X_weights %*% (
    post_chain$post_eta[1:dim_Xw, ]
  )),
  t(X_weights %*% (
    post_chain$post_eta[(dim_Xw + 1):(dim_Xw * 2), ]
  ))))
  for (g in 3:(n_group - 1)) {
    c_eta_X = cbind(c_eta_X, t(pnorm(X_weights %*% (
      post_chain$post_eta[((g - 1) * dim_Xw + 1):(g * dim_Xw), ]
    ))))
  }
  c_eta_X = cbind(c_eta_X, 1)
  
  t = cbind(c_eta_X[, 1],
            c_eta_X[, 2] * (1 - c_eta_X[, 1]),
            sapply(3:n_group, function(g)
              c_eta_X[, g] * apply(1 - c_eta_X[, 1:(g - 1)], 1, prod)))
  
  clusters = apply(t, 1, which.max)
  
  beta = post_chain$post_alpha[seq(2, dim(post_chain$post_alpha)[1], dim_REG), ]
  intercept = post_chain$post_alpha[-seq(2, dim(post_chain$post_alpha)[1], dim_REG), ]
  
  if (dim_REG < 3) {
    values = sapply(1:length(clusters), function(c)
      intercept[clusters[c], c] + x %*% beta[clusters[c], c])
  } else{
    values = sapply(1:length(clusters), function(c)
      t(apply(X_tilde[, -2], 2, mean)) %*% intercept[((clusters[c] - 1) *
                                                        (dim_REG - 1) + 1):(clusters[c] * (dim_REG - 1)), c] +
        x %*% beta[clusters[c], c])
  }
  
  #return(quantile(values, prob=0.5, na.rm=TRUE))
  return(values)
}

curve_ADJ <- function(x, post_chain, X_tilde) {
  X_weights = split_4quantile_x(data_analysis$X, x)
  dim_Xw = dim(X_weights)[2]
  dim_REG = dim(X_tilde)[2]
  
  c_eta_X = pnorm(cbind(t(X_weights %*% (
    post_chain$post_eta[1:dim_Xw, ]
  )),
  t(X_weights %*% (
    post_chain$post_eta[(dim_Xw + 1):(dim_Xw * 2), ]
  ))))
  for (g in 3:(n_group - 1)) {
    c_eta_X = cbind(c_eta_X, t(pnorm(X_weights %*% (
      post_chain$post_eta[((g - 1) * dim_Xw + 1):(g * dim_Xw), ]
    ))))
  }
  c_eta_X = cbind(c_eta_X, 1)
  
  t = cbind(c_eta_X[, 1],
            c_eta_X[, 2] * (1 - c_eta_X[, 1]),
            sapply(3:n_group, function(g)
              c_eta_X[, g] * apply(1 - c_eta_X[, 1:(g - 1)], 1, prod)))
  
  clusters = apply(t, 1, which.max)
  
  values = sapply(1:length(clusters), function(c)
    post_chain$post_intercept[clusters[c], c] + x %*% post_chain$post_beta_x[clusters[c], c])
  
  #return(quantile(values, prob=0.5, na.rm=TRUE))
  return(values)
}


##########  estimation curves ###############

min_x = min(data_analysis$X)
min_x = 5
max_x = max(data_analysis$X)
max_x = 14

# causal effect for a grid of points
points_x_i = seq(min_x, max_x, length.out = 100)

points_YX = sapply(points_x_i, function(x)
  curve_ADJ(x, post_chain = model_YX, X_tilde = X_tilde_YX))
points_YXC = sapply(points_x_i, function(x)
  curve_ADJ(x, post_chain = model_YXC, X_tilde = X_tilde_YXC))
points_YXU = sapply(points_x_i, function(x)
  curve_ADJ(x, post_chain = model_YXU, X_tilde = X_tilde_YXU))
points_YXCU = sapply(points_x_i, function(x)
  curve_ADJ(x, post_chain = model_YXCU, X_tilde = X_tilde_YXCU))
points_YXCWZ = sapply(points_x_i, function(x)
  curve_ADJ(x, post_chain = model_YXCWZ, X_tilde = X_tilde_YXCWZ))

points_adj = sapply(points_x_i, function(x)
  curve_ADJ(x, post_chain = model_adj, X_tilde = X_tilde_adj))
points_adjC = sapply(points_x_i, function(x)
  curve_ADJ(x, post_chain = model_adjC, X_tilde = X_tilde_adjC))
points_adjCW = sapply(points_x_i, function(x)
  curve_ADJ(x, post_chain = model_adjCW, X_tilde = X_tilde_adjC))

smooth_YX = ksmooth(points_x_i,
                    apply(points_YX, 2, quantile, prob = 0.5, na.rm = TRUE),
                    "normal",
                    bandwidth = 0.5)$y
smooth_YXC = ksmooth(
  points_x_i,
  apply(points_YXC, 2, quantile, prob = 0.5, na.rm = TRUE),
  "normal",
  bandwidth = 0.5
)$y
smooth_YXU = ksmooth(
  points_x_i,
  apply(points_YXU[1:600, ], 2, quantile, prob = 0.5, na.rm =
          TRUE),
  "normal",
  bandwidth = 0.5
)$y
smooth_YXCU = ksmooth(
  points_x_i,
  apply(points_YXCU, 2, quantile, prob = 0.5, na.rm =
          TRUE),
  "normal",
  bandwidth = 0.5
)$y
smooth_YXCWZ = ksmooth(
  points_x_i,
  apply(points_YXCWZ, 2, quantile, prob = 0.5, na.rm = TRUE),
  "normal",
  bandwidth = 0.5
)$y

smooth_YX = ksmooth(points_x_i,
                    apply(points_YX, 2, quantile, prob = 0.5, na.rm = TRUE),
                    "normal",
                    bandwidth = 0.5)$y
smooth_adjC = ksmooth(
  points_x_i,
  apply(points_adjC, 2, quantile, prob = 0.5, na.rm =
          TRUE),
  "normal",
  bandwidth = 0.5
)$y
smooth_adjCW = ksmooth(
  points_x_i,
  apply(points_adjCW, 2, quantile, prob = 0.5, na.rm =
          TRUE),
  "normal",
  bandwidth = 0.5
)$y

save.image("results_models_data_2.RData")
