# 4 splits : used in simulations section
split_4quantile <- function(x_matrix) {
  # splitting points
  seq_X = c(min(x_matrix), quantile(x_matrix, c(0.25, 0.5, 0.75)), max(x_matrix))
  
  # creating matrix for the probit regression in the weights
  X_weights = cbind(
    cbind(1, x_matrix) * I(x_matrix < seq_X[2]),
    cbind(1, x_matrix) * I(x_matrix < seq_X[3] &
                             x_matrix > seq_X[2]),
    cbind(1, x_matrix) * I(x_matrix < seq_X[4] &
                             x_matrix > seq_X[3]),
    cbind(1, x_matrix) * I(x_matrix > seq_X[4])
  )
  X_weights = (X_weights != 0) * (X_weights) +
    (X_weights == 0) * matrix(rnorm(prod(dim(X_weights)), 0, 0.003),
                              ncol = dim(X_weights)[2],
                              nrow = dim(X_weights)[1])
  
  return(X_weights)
  
}

# estimating CERF with our proposed methods
curve_ADJ <- function(x,
                      post_chain,
                      X_tilde,
                      x_split_method,
                      probs = NULL,
                      pts = NULL,
                      data_analysis,
                      n_group) {
  # preparing matrix for the probit regression in the weights
  if (x_split_method == "split_quantile_x") {
    X_weights = split_quantile_x(data_analysis$X, x, probs = probs)
  } else if (x_split_method == "split_at_fixed_pt_x") {
    X_weights = split_at_fixed_pt_x(data_analysis$X, x, pts = pts)
  } else{
    X_weights = split_4quantile_x(data_analysis$X, x)
  }
  
  # dimensions matrices
  dim_Xw = dim(X_weights)[2]
  dim_REG = dim(X_tilde)[2]
  
  # estimation weights for a grid of values of treatment X
  c_eta_X = pnorm(cbind(t(X_weights %*% (
    post_chain$post_eta[1:dim_Xw,]
  )),
  t(X_weights %*% (
    post_chain$post_eta[(dim_Xw + 1):(dim_Xw * 2),]
  ))))
  for (g in 3:(n_group - 1)) {
    c_eta_X = cbind(c_eta_X, t(pnorm(X_weights %*% (
      post_chain$post_eta[((g - 1) * dim_Xw + 1):(g * dim_Xw),]
    ))))
  }
  c_eta_X = cbind(c_eta_X, 1)
  t = cbind(c_eta_X[, 1],
            c_eta_X[, 2] * (1 - c_eta_X[, 1]),
            sapply(3:n_group, function(g)
              c_eta_X[, g] * apply(1 - c_eta_X[, 1:(g - 1)], 1, prod)))
  # component mixture allocation
  clusters = apply(t, 1, which.max)
  # estimation CERF for the grid of values for X
  values = sapply(1:length(clusters), function(c)
    post_chain$post_intercept[clusters[c], c] + x %*% post_chain$post_beta_x[clusters[c], c])
  
  #alternative : return(quantile(values, prob=0.5, na.rm=TRUE))
  return(values)
}
