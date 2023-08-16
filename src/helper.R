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



split_4quantile2 <- function(x_matrix) {
  seq_X = c(min(x_matrix),7,  10, 13, max(x_matrix))
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


split_4quantile2_x <- function(x_matrix, x) {
  seq_X = c(min(x_matrix), 7, 10, 13, max(x_matrix))
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

split_5quantile <- function(x_matrix) {
  seq_X = c(min(x_matrix), quantile(x_matrix, c(0.1, 0.25, 0.5, 0.75)), max(x_matrix))
  X_weights = cbind(
    x_matrix * I(x_matrix < seq_X[2]),
    x_matrix * I(x_matrix < seq_X[3] &
                   x_matrix > seq_X[2]),
    x_matrix * I(x_matrix < seq_X[4] &
                   x_matrix > seq_X[3]),
    x_matrix * I(x_matrix < seq_X[5] &
                   x_matrix > seq_X[4]),
    x_matrix * I(x_matrix > seq_X[5])
  )
  X_weights = (X_weights != 0) * (X_weights) +
    (X_weights == 0) * matrix(rnorm(prod(dim(X_weights)), 0, 0.003),
                              ncol = dim(X_weights)[2],
                              nrow = dim(X_weights)[1])
  
  return(X_weights)
}


split_5quantile_x <- function(x_matrix, x) {
  seq_X = c(min(x_matrix), quantile(x_matrix, c(0.1, 0.25, 0.5, 0.75)), max(x_matrix))
  X_weights = cbind(
    x * I(x < seq_X[2]),
    x * I(x < seq_X[3] & x > seq_X[2]),
    x * I(x < seq_X[4] & x > seq_X[3]),
    x * I(x < seq_X[5] & x > seq_X[4]),
    x * I(x > seq_X[5])
  )
  X_weights = (X_weights != 0) * (X_weights) +
    (X_weights == 0) * matrix(rnorm(prod(dim(X_weights)), 0, 0.003),
                              ncol = dim(X_weights)[2],
                              nrow = dim(X_weights)[1])
  
  return(X_weights)
}


split_5quantile2 <- function(x_matrix) {
  seq_X = c(min(x_matrix), quantile(x_matrix, c(0.2, 0.4, 0.6, 0.8)), max(x_matrix))
  X_weights = cbind(
    x_matrix * I(x_matrix < seq_X[2]),
    x_matrix * I(x_matrix < seq_X[3] &
                   x_matrix > seq_X[2]),
    x_matrix * I(x_matrix < seq_X[4] &
                   x_matrix > seq_X[3]),
    x_matrix * I(x_matrix < seq_X[5] &
                   x_matrix > seq_X[4]),
    x_matrix * I(x_matrix > seq_X[5])
  )
  X_weights = (X_weights != 0) * (X_weights) +
    (X_weights == 0) * matrix(rnorm(prod(dim(X_weights)), 0, 0.003),
                              ncol = dim(X_weights)[2],
                              nrow = dim(X_weights)[1])
  
  return(X_weights)
}


split_5quantile2_x <- function(x_matrix, x) {
  seq_X = c(min(x_matrix), quantile(x_matrix, c(0.2, 0.4, 0.6, 0.8)), max(x_matrix))
  X_weights = cbind(
    x * I(x < seq_X[2]),
    x * I(x < seq_X[3] & x > seq_X[2]),
    x * I(x < seq_X[4] & x > seq_X[3]),
    x * I(x < seq_X[5] & x > seq_X[4]),
    x * I(x > seq_X[5])
  )
  X_weights = (X_weights != 0) * (X_weights) +
    (X_weights == 0) * matrix(rnorm(prod(dim(X_weights)), 0, 0.003),
                              ncol = dim(X_weights)[2],
                              nrow = dim(X_weights)[1])
  
  return(X_weights)
}


split_6quantile <- function(x_matrix) {
  seq_X = c(min(x_matrix), quantile(x_matrix, c(0.15, 0.3, 0.45, 0.6, 0.75)), max(x_matrix))
  X_weights = cbind(
    x_matrix * I(x_matrix < seq_X[2]),
    x_matrix * I(x_matrix < seq_X[3] &
                   x_matrix > seq_X[2]),
    x_matrix * I(x_matrix < seq_X[4] &
                   x_matrix > seq_X[3]),
    x_matrix * I(x_matrix < seq_X[5] &
                   x_matrix > seq_X[4]),
    x_matrix * I(x_matrix < seq_X[6] &
                   x_matrix > seq_X[5]),
    x_matrix * I(x_matrix > seq_X[6])
  )
  X_weights = (X_weights != 0) * (X_weights) +
    (X_weights == 0) * matrix(rnorm(prod(dim(X_weights)), 0, 0.003),
                              ncol = dim(X_weights)[2],
                              nrow = dim(X_weights)[1])
  
  return(X_weights)
}


split_6quantile_x <- function(x_matrix, x) {
  seq_X = c(min(x_matrix), quantile(x_matrix, c(0.15, 0.3, 0.45, 0.6, 0.75)), max(x_matrix))
  X_weights = cbind(
    x * I(x < seq_X[2]),
    x * I(x < seq_X[3] & x > seq_X[2]),
    x * I(x < seq_X[4] & x > seq_X[3]),
    x * I(x < seq_X[5] & x > seq_X[4]),
    x * I(x < seq_X[6] & x > seq_X[6]),
    x * I(x > seq_X[6])
  )
  X_weights = (X_weights != 0) * (X_weights) +
    (X_weights == 0) * matrix(rnorm(prod(dim(X_weights)), 0, 0.003),
                              ncol = dim(X_weights)[2],
                              nrow = dim(X_weights)[1])
  
  return(X_weights)
}

