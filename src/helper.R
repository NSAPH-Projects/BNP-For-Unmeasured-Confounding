###################################################################################
# ---   Tools to Increase Flexibility of Bayesian Nonparametric Models           ---
###################################################################################

# These functions divide exposure/treatment X to different segments, create
# design matrix. By this means, at the weights construction steps, we allow
# coefficients in linear regression  models to be different for different
# segments of X

# split training data X based on quantile
split_4quantile <- function(x_matrix) {
  seq_X = c(min(x_matrix), quantile(x_matrix, c(0.25, 0.5, 0.75)), max(x_matrix))
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

# split testing data x based on quantile
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




# split training data X based on a sequence of given points
split_at_fixed_pt <- function(x_matrix, pts = c(7, 10, 13)) {
  n = length(x_matrix)
  seq_X = c(min(x_matrix), pts, max(x_matrix))
  X_weights = matrix(0, nrow = n, ncol = length(seq_X) - 1)
  for (j in 1:n) {
    for (i in 1:(length(seq_X) - 1)) {
      if (x_matrix[j] >= seq_X [i] && x_matrix[j] < seq_X [i + 1]) {
        X_weights[j, i] = x_matrix[j]
      }
    }
  }
  
  X_weights = (X_weights != 0) * (X_weights) +
    (X_weights == 0) * matrix(rnorm(prod(dim(X_weights)), 0, 0.003),
                              ncol = dim(X_weights)[2],
                              nrow = dim(X_weights)[1])
  
  return(X_weights)
}

# split testing data x based on a sequence of given points
split_at_fixed_pt_x <- function(x_matrix, x, pts = c(7, 10, 13)) {
  seq_X = c(min(x_matrix), pts, max(x_matrix))
  n = length(x)
  X_weights = matrix(0, nrow = n, ncol = length(seq_X) - 1)
  seq_X = c(min(x_matrix), pts, max(x_matrix))
  for (j in 1:n) {
    for (i in 1:(length(seq_X) - 1)) {
      if (x[j] >= seq_X [i] && x[j] < seq_X [i + 1]) {
        X_weights[j, i] = x[j]
      }
    }
  }
  X_weights = (X_weights != 0) * (X_weights) +
    (X_weights == 0) * matrix(rnorm(prod(dim(X_weights)), 0, 0.003),
                              ncol = dim(X_weights)[2],
                              nrow = dim(X_weights)[1])
  
  return(X_weights)
}

#  split training data X for quantiles corresponding to the given probabilities
split_quantile <- function(x_matrix, probs = c(0.2, 0.4, 0.6, 0.8)) {
  # example: split_quantile(x_matrix, probs = c(0.1, 0.25, 0.5, 0.75))
  #        : split_quantile(x_matrix, probs = c(0.15, 0.3, 0.45, 0.6, 0.75))
  seq_X = c(min(x_matrix),
            quantile(x_matrix, probs = probs),
            max(x_matrix))
  n = length(x_matrix)
  X_weights = matrix(0, nrow = n, ncol = length(seq_X) - 1)
  seq_X = c(min(x_matrix), probs, max(x_matrix))
  for (j in 1:n) {
    for (i in 1:(length(seq_X) - 1)) {
      if (x_matrix[j] >= seq_X [i] && x_matrix[j] < seq_X [i + 1]) {
        X_weights[j, i] = x_matrix[j]
      }
    }
  }
  X_weights = (X_weights != 0) * (X_weights) +
    (X_weights == 0) * matrix(rnorm(prod(dim(X_weights)), 0, 0.003),
                              ncol = dim(X_weights)[2],
                              nrow = dim(X_weights)[1])
  
  return(X_weights)
}

#  split testing data x for quantiles corresponding to the given probabilities
split_quantile_x <-
  function(x_matrix, x, probs = c(0.2, 0.4, 0.6, 0.8)) {
    # example: split_quantile_x(x_matrix, x, probs = c(0.1, 0.25, 0.5, 0.75))
    #        : split_quantile_x(x_matrix, x, probs = c(0.15, 0.3, 0.45, 0.6, 0.75))
    seq_X = c(min(x_matrix),
              quantile(x_matrix, probs = probs),
              max(x_matrix))
    n = length(x)
    X_weights = matrix(0, nrow = n, ncol = length(seq_X) - 1)
    seq_X = c(min(x_matrix), probs, max(x_matrix))
    for (j in 1:n) {
      for (i in 1:(length(seq_X) - 1)) {
        if (x[j] >= seq_X [i] && x[j] < seq_X [i + 1]) {
          X_weights[j, i] = x[j]
        }
      }
    }
    X_weights = (X_weights != 0) * (X_weights) +
      (X_weights == 0) * matrix(rnorm(prod(dim(X_weights)), 0, 0.003),
                                ncol = dim(X_weights)[2],
                                nrow = dim(X_weights)[1])
    
    return(X_weights)
  }
