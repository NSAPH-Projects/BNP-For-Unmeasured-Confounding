# gibbs sampler used in simulations section
DDP_ADJ <- function(s_seed = 1,
                    X_tilde,
                    X_W_reg,
                    X_split_method,
                    probs = NULL,
                    pts = NULL,
                    R = R,
                    R_burnin,
                    n_group,
                    n = n,
                    data_analysis) {
  # seed for reproducibility
  set.seed(s_seed)
  
  # ------   preparing variables   ------
  
  # matrix for probit regressio in the weights
  dim_REG = dim(X_tilde)[2]
  dim_W = dim(X_W_reg)[2]
  if (X_split_method == "split_quantile") {
    X_weights = split_quantile(data_analysis$X, probs = probs)
  } else if (X_split_method == "split_at_fixed_pt") {
    X_weights = split_at_fixed_pt(data_analysis$X, pts = pts)
  } else if (X_split_method == "split_4quantile") {
    X_weights = split_4quantile(data_analysis$X)
  } else if (X_split_method == "split_4quantile_centered") {
    X_weights = split_4quantile_centered(data_analysis$X)
  } else {
    print ("not supported")
  }
  
  # ------   hyperparameters   -----
  
  # alpha = parameters of the regression mean in the components of Y|X,Z mixture
  # mu(=mean) and sigma(=var) for the normal prior for alpha
  alpha_mu = rep(0, dim_REG)
  alpha_sigma = 1
  
  # delta = parameters of the mean in the distribution of W|X,Z
  delta_mu = rep(0, dim_W)
  delta_sigma = 1
  
  # parameters for gamma_y \sim InvGamma
  # gamma_y = variance of the components of Y|X,Z mixture
  gamma_y_prior = c(1, 0.5)
  
  # parameters for gamma_w \sim InvGamma
  # gamma_w = variance of the distribution of W|X,Z
  gamma_w_prior = c(2, 2)
  # eta's iperparameters
  # eta = parameters in the probit regression in the weights
  eta_prior = c(0, 5)
  
  # ------   initialitation   -----
  
  # parameters
  alpha = matrix(rep(rmvnorm(
    1, alpha_mu * rep(1, dim_REG),
    alpha_sigma * diag(dim_REG)
  ), n_group),
  ncol = n_group, nrow = dim_REG)
  delta = rnorm(dim_W, delta_mu, delta_sigma)
  gamma_y = rep(rinvgamma(1, gamma_y_prior[1], gamma_y_prior[2]), n_group)
  gamma_w = rinvgamma(1, gamma_w_prior[1], gamma_w_prior[2])
  eta = matrix(
    rnorm((n_group - 1) *  dim(X_weights)[2], eta_prior[1], eta_prior[2]),
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
  
  # cluster allocation variables
  K = sample(1:n_group, n, replace = TRUE)
  table_Ki = rep(0, n_group)
  table_Ki[sort(unique(K))] = table(K)
  # latent variable for augmentation scheme
  Q = matrix(NA, nrow = n, ncol = n_group - 1)
  # causal effect parameter
  beta_x = rep(0, n_group)
  
  # ------   saving informations   -----
  
  # empty matrix where save all the informations for each iteration
  post_alpha = matrix(NA, nrow = n_group * dim_REG, ncol = R)
  post_delta = matrix(NA, nrow = dim_W, ncol = R)
  post_beta_x = matrix(NA, nrow = n_group, ncol = R)
  post_intercept = matrix(NA, nrow = n_group, ncol = R)
  post_eta = matrix(NA, nrow = (n_group - 1) * (dim(X_weights)[2]), ncol =
                      R)
  K_all = matrix(NA, nrow = n, ncol = R)
  
  for (r in 1:R) {
    #######################################################
    # -----  estimation parameter for W|X,Z model  ------
    #######################################################
    
    #delta  ( parameters of the mean in the distribution of W|X,Z )
    V_delta = diag(dim_W) / delta_sigma + t(X_W_reg) %*% X_W_reg / gamma_w
    m_delta = delta_mu / delta_sigma + t(X_W_reg) %*% data_analysis$W /
      gamma_w
    delta = rmvnorm(1, solve(V_delta) %*% m_delta, solve(V_delta))
    
    #gamma_w    ( parameter of the variance in the distribution of W|X,Z )
    g_w = sum((data_analysis$W - X_W_reg %*% t(delta)) ^ 2)
    gamma_w = rinvgamma(1, gamma_w_prior[1] + n / 2, gamma_w_prior[2] +
                          g_w / 2)
    
    #######################################################
    # -----  estimation parameter for Y|X,Z model  ------
    #######################################################
    # ---- 1: Mixture Component Specific Parameters  ----
    #######################################################
    
    #for compute eta (regression parameters in the weights)
    # we need to compute mu_Q and Q
    
    # mu_Q: mean for the gaussian variable Q
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
    
    # augmentation data: variable Q
    for (i in 1:n) {
      for (l in 1:(min(K[i], n_group - 1))) {
        if (l < K[i]) {
          Q[i, l] = rtruncnorm(1, b = 0, mean = qnorm(mu_Q[i, l]))
        } else{
          Q[i, l] = rtruncnorm(1, a = 0, mean = qnorm(mu_Q[i, l]))
        }
      }
    }
    
    # eta = parameters in the probit regression in the weights
    check_q = sapply(1:(n_group - 1), function(l)
      length(Q[K >= l, l]))
    v_eta = lapply(which(check_q > 3), function(l)
      1 / eta_prior[2] + t(X_weights[K >= l,]) %*% X_weights[K >= l,])
    m_eta = sapply(which(check_q > 3), function(l)
      eta_prior[1] / eta_prior[2] +
        t(X_weights[K >= l,]) %*% Q[K >= l, l] / gamma_y[l])
    # each component has to have at least 3 units to estimate eta
    if (any(check_q < 4)) {
      for (l in which(check_q == 1 | check_q == 2 | check_q == 3)) {
        v_eta[[l]] = diag(dim(X_weights)[2]) / eta_prior[2]
        m_eta = cbind(m_eta, rep(eta_prior[1] / eta_prior[2], dim(X_weights)[2]))
      }
    }
    eta[, which(table_Ki[-n_group] > 0)] = sapply(which(table_Ki[-n_group] >
                                                          0),
                                                  function(l)
                                                    rmvnorm(1, solve(v_eta[[l]]) %*% m_eta[, l]))
    eta[, which(table_Ki[-n_group] == 0)] = rnorm(length(which(table_Ki[-n_group] ==
                                                                 0)) * dim(X_weights)[2],
                                                  eta_prior[1],
                                                  sqrt(eta_prior[2]))
    
    
    ##############################################
    # -------- 2: Component Allocation  ---------
    ##############################################
    
    #  tau = recursive weigths in the mixture
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
      tau[i,] = 1 / n_group
    }
    
    # K = allocation component variables
    # it says each unit in which component of the mixture is allocated
    K = sapply(1:n, function(i)
      (1:n_group) %*% rmultinom(1, 1, tau[i,]))
    table_Ki = rep(0, n_group)
    table_Ki[sort(unique(K))] = table(K)
    
    ##############################################
    # ------- 3: Regression mean of Y|X,Z  ------
    ##############################################
    
    # alpha + gamma_y :
    # parameters of the mean + variance in each mixture components of the distribution of Y|X,Z
    for (g_k in 1:n_group) {
      if (length(which(K == g_k)) == 1) {
        V = diag(dim_REG) / alpha_sigma + t(t(X_tilde[K == g_k,])) %*% t(X_tilde[K ==
                                                                                   g_k,]) / gamma_y[g_k]
        M = alpha_mu / alpha_sigma + t(t(X_tilde[K == g_k,])) %*% t(data_analysis$Y[K ==
                                                                                      g_k]) / gamma_y[g_k]
      } else{
        V = diag(dim_REG) / alpha_sigma + t(X_tilde[K == g_k,]) %*% X_tilde[K ==
                                                                              g_k,] / gamma_y[g_k]
        M = alpha_mu / alpha_sigma + t(X_tilde[K == g_k,]) %*% data_analysis$Y[K ==
                                                                                 g_k] / gamma_y[g_k]
      }
      alpha[, g_k] = t(rmvnorm(1, solve(V) %*% M, solve(V)))
      
      G = sum((data_analysis$Y[K == g_k] - X_tilde[K == g_k,] %*% (alpha[, g_k])) ^
                2)
      gamma_y[g_k] = rinvgamma(1, gamma_y_prior[1] + (sum(K == g_k)) / 2, gamma_y_prior[2] +
                                 G / 2)
    }
    
    
    ##############################################
    # --- Causal effects: parameters for CERF  ----
    ##############################################
    
    #beta_x + intercept
    beta_x = alpha[2,] - alpha[3,] * delta[1, 2] / delta[1, 3]
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
    
    # -----   saving information   -----
    # parameters
    post_alpha[, r] = c(alpha)
    post_delta[, r] = delta[, 1:dim_W]
    post_eta[, r] = c(eta)
    # CERF parameters
    post_beta_x[, r] = beta_x
    post_intercept[, r] = intercept
    
    #K membership
    K_all[, r] = K
    
    #
    if (r %% 500 == 0)
      print(paste0(r, "/", R, " iterations"))
  }
  
  return(
    list(
      post_alpha = post_alpha[, (R_burnin + 1):R],
      post_beta_x = post_beta_x[, (R_burnin + 1):R],
      post_intercept = post_intercept[, (R_burnin + 1):R],
      post_eta = post_eta[, (R_burnin + 1):R],
      post_delta = post_delta[, (R_burnin + 1):R],
      K_all = K_all[, (R_burnin + 1):R]
    )
  )
}