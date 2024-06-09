# 4 splits : used in simulations section
split_4quantile <- function(x_matrix) {
  # splitting points
  seq_X = c(min(x_matrix), quantile(x_matrix, c(0.25, 0.5, 0.75)), max(x_matrix))
  
  # creating matrix for the probit regression in the weights
  X_weights = cbind(
    x_matrix * I(x_matrix < seq_X[2]),
    x_matrix * I(x_matrix < seq_X[3] &
                   x_matrix > seq_X[2]),
    x_matrix * I(x_matrix < seq_X[4] &
                   x_matrix > seq_X[3]),
    x_matrix * I(x_matrix > seq_X[4])
  )
  X_weights = (X_weights != 0) * (X_weights) +
    (X_weights == 0) * matrix(rnorm(prod(dim(X_weights)), 0, 0.003),
                              ncol = dim(X_weights)[2],
                              nrow = dim(X_weights)[1])
  
  return(X_weights)
}

# estimating CERF with our proposed methods
curve_chains_with <-
  function(x, data_sim, post_chian, split_function_x) {
    # matrix with the point in the support of X (=treatment) considered
    X_weights = split_function_x(data_sim$X, x)
    
    # expected weights for the points x considered
    c_eta_X = pnorm(cbind(t(X_weights %*% (
      post_chian$post_eta[1:(dim(X_weights)[2]), ]
    )),
    t(X_weights %*% (
      post_chian$post_eta[(dim(X_weights)[2] + 1):(dim(X_weights)[2] * 2), ]
    ))))
    for (g in 3:(n_group - 1)) {
      c_eta_X = cbind(c_eta_X, t(pnorm(X_weights %*% (
        post_chian$post_eta[((g - 1) * (dim(X_weights)[2]) + 1):(g * (dim(X_weights)[2])), ]
      ))))
    }
    c_eta_X = cbind(c_eta_X, 1)
    
    # expected value of CERF for each value of x
    t = cbind(c_eta_X[, 1],
              c_eta_X[, 2] * (1 - c_eta_X[, 1]),
              sapply(3:n_group, function(g)
                c_eta_X[, g] * apply(1 - c_eta_X[, 1:(g - 1)], 1, prod)))
    clusters = apply(t, 1, which.max)
    values = sapply(1:length(clusters), function(c)
      post_chian$post_intercept[clusters[c], c] +
        x %*% post_chian$post_beta_x[clusters[c], c])
    
    return(quantile(values, prob = 0.5, na.rm = TRUE))
  }
