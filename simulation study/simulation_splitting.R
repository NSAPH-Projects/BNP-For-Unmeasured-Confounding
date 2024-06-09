###################################################################################
# ---              SIMULATION STUDY             ---
# ---                4 settings                 ---
# ---              5000 sample size             ---
###################################################################################

# library
library(mvtnorm)
library(invgamma)

###################################################################################
###################################################################################
###################################################################################
# continuous functions for Y (outcome)
# given treatment (X) and U (unmeasured confounder)

# --- MIXTURE OF LINEARS ---
fun_Y_1 <- function(X, U) {
  return(I(X < 5.5) * (+1 + 2 * X + 2 * U) + I(X >= 5.5) * (-16 + 5 * X +
                                                              2.5 * U))
}
fun_Y_1_int <- function(X, U) {
  return(fun_Y_1(X, mean(U)))
}

# --- PARABOLA ---
fun_Y_2 <- function(X, U) {
  return(-10 + 1.5 * (X - 6) ^ 2 + 4 * U)
}
fun_Y_2_int <- function(X, U) {
  return(fun_Y_2(X, mean(U)))
}

# --- S FUNCTION ---
fun_Y_3 <- function(X, U) {
  return(1 / (1 + exp(-5 * (X - 5))) + 1.7 * U)
}
fun_Y_3_int <- function(X, U) {
  return(fun_Y_3(X, mean(U)))
}

# --- HALF HORIZONTAL PARABOLA ---
fun_Y_4 <- function(X, U) {
  return(-2 * exp(-1.4 * (X - 6)) + 3 * exp(U))
}
fun_Y_4_int <- function(X, U) {
  return(-2 * exp(-1.4 * (X - 6)) + 3 * mean(exp((U))))
}
###################################################################################
# general informations

n = 6250       #number of units for each sample, corespond to 5000 points
#(I will cut the tails because too scattered)

###################################################################################

# --- generation of samplers for the 4 settings ---

# generating functions
simulation_data_1 <- function(n) {
  U = rnorm(n, 1, 0.2)
  Z = rnorm(n,-1 + 1.5 * U, 0.2)
  W = rnorm(n, 1 - 2 * U, 0.2)
  X = 1.5 + 4 * U + rnorm(n, 0, 0.2)
  Y = rnorm(n, fun_Y_1(X, U), 0.3)
  a = cbind(U, W, Z, X, Y)
  a = a[X > quantile(X, 0.1) & X < quantile(X, 0.9), ]
  return(as.data.frame(a))
}

simulation_data_2 <- function(n) {
  U = rnorm(n, 1, 0.2)
  Z = rnorm(n,-1 + 1.5 * U, 0.2)
  W = rnorm(n, 1 - 2 * U, 0.2)
  X = 2.5 + 4 * U + rnorm(n, 0, 0.2)
  Y = rnorm(n, fun_Y_2(X, U), 0.2)
  a = cbind(U, W, Z, X, Y)
  a = a[X > quantile(X, 0.1) & X < quantile(X, 0.9), ]
  return(as.data.frame(a))
}

simulation_data_3 <- function(n) {
  U = rnorm(n, 1, 0.2)
  Z = rnorm(n,-1 + 1.5 * U, 0.2)
  W = rnorm(n, 1 - 2 * U, 0.2)
  X = 1 + 4 * U + rnorm(n, 0, 0.2)
  Y = rnorm(n, fun_Y_3(X, U), 0.1)
  a = cbind(U, W, Z, X, Y)
  a = a[X > quantile(X, 0.1) & X < quantile(X, 0.9), ]
  return(as.data.frame(a))
}

simulation_data_4 <- function(n) {
  U = rnorm(n, 1, 0.2)
  Z = rnorm(n,-1 + 1.5 * U, 0.2)
  W = rnorm(n, 1 - 2 * U, 0.2)
  X = 2.5 + 4 * U + rnorm(n, 0,  0.2)
  Y = rnorm(n, fun_Y_4(X, U), 0.2)
  a = cbind(U, W, Z, X, Y)
  a = a[X > quantile(X, 0.1) & X < quantile(X, 0.9), ]
  return(as.data.frame(a))
}

# --- Simulation ---
n = 6250
data_sim_1 = simulation_data_1(n)
n = dim(data_sim_1)[1]
X_tilde_adj = as.matrix(cbind(rep(1, n), data_sim_1$X, data_sim_1$Z))
X_split_method = "split_4quantile"

model_adj1 <- DDP_ADJ(
  s_seed = 10,
  X_tilde = X_tilde_adj,
  X_W_reg = X_tilde_adj,
  X_split_method,
  probs = NULL,
  pts = NULL,
  R = 4000,
  R_burnin = 3000,
  n_group = 30,
  n = n,
  data_sim_1
)

plot_posterior_splitting(
  data_sim_1,
  model_adj1,
  X_tilde = X_tilde_adj,
  x_split_method = "split_4quantile_x",
  probs = NULL,
  pts = NULL,
  n_group = 30,
  fun_Y_1_int,
  bandwidth = 0.2,
  smooth = TRUE
)

pdf(
  "splitting_smoothing_n=5000_setting_A_senario1.pdf",
  width = 9,
  height = 7
)
ggplot_posterior_splitting(
  data_sim_1,
  model_adj1,
  X_tilde = X_tilde_adj,
  x_split_method = "split_4quantile_x",
  probs = NULL,
  pts = NULL,
  n_group = 30,
  funct_YU = fun_Y_1_int,
  bandwidth = 0.2,
  smooth_curve = TRUE,
  num = 1
)
dev.off()

pdf("splitting_n=5000_setting_A_senario1.pdf",
    width = 9,
    height = 7)
ggplot_posterior_splitting(
  data_sim_1,
  model_adj1,
  X_tilde = X_tilde_adj,
  x_split_method = "split_4quantile_x",
  probs = NULL,
  pts = NULL,
  n_group = 30,
  funct_YU = fun_Y_1_int,
  bandwidth = 0.2,
  smooth_curve = FALSE,
  num = 1
)
dev.off()


n = 6250
data_sim_2 = simulation_data_2(n)
n = dim(data_sim_2)[1]
X_tilde_adj = as.matrix(cbind(rep(1, n), data_sim_2$X, data_sim_2$Z))
X_split_method = "split_4quantile"

model_adj2 <- DDP_ADJ(
  s_seed = 10,
  X_tilde = X_tilde_adj,
  X_W_reg = X_tilde_adj,
  X_split_method,
  probs = NULL,
  pts = NULL,
  R = 4000,
  R_burnin = 3000,
  n_group = 30,
  n = n,
  data_sim_2
)

plot_posterior_splitting(
  data_sim_2,
  model_adj,
  X_tilde = X_tilde_adj,
  x_split_method = "split_4quantile_x",
  probs = NULL,
  pts = NULL,
  n_group = 30,
  fun_Y_2_int
)

pdf(
  "splitting_smoothing_n=5000_setting_A_senario2.pdf",
  width = 9,
  height = 7
)
ggplot_posterior_splitting(
  data_sim_2,
  model_adj2,
  X_tilde = X_tilde_adj,
  x_split_method = "split_4quantile_x",
  probs = NULL,
  pts = NULL,
  n_group = 30,
  funct_YU = fun_Y_2_int,
  bandwidth = 0.2,
  smooth_curve = TRUE,
  num = 2
)
dev.off()

pdf("splitting_n=5000_setting_A_senario2.pdf",
    width = 9,
    height = 7)
ggplot_posterior_splitting(
  data_sim_2,
  model_adj2,
  X_tilde = X_tilde_adj,
  x_split_method = "split_4quantile_x",
  probs = NULL,
  pts = NULL,
  n_group = 30,
  funct_YU = fun_Y_2_int,
  bandwidth = 0.2,
  smooth_curve = FALSE,
  num = 2
)
dev.off()


n = 6250
data_sim_3 = simulation_data_3(n)
n = dim(data_sim_3)[1]
X_tilde_adj = as.matrix(cbind(rep(1, n), data_sim_3$X, data_sim_3$Z))
X_split_method = "split_4quantile"

model_adj3 <- DDP_ADJ(
  s_seed = 5,
  X_tilde = X_tilde_adj,
  X_W_reg = X_tilde_adj,
  X_split_method,
  probs = NULL,
  pts = NULL,
  R = 4000,
  R_burnin = 3000,
  n_group = 30,
  n = n,
  data_sim_3
)

plot_posterior_splitting(
  data_sim_3,
  model_adj3,
  X_tilde = X_tilde_adj,
  x_split_method = "split_4quantile_x",
  probs = NULL,
  pts = NULL,
  n_group = 30,
  fun_Y_3_int,
  bandwidth = 0.2,
  smooth = TRUE
)

pdf(
  "splitting_smoothing_n=5000_setting_A_senario3.pdf",
  width = 9,
  height = 7
)
ggplot_posterior_splitting(
  data_sim_3,
  model_adj3,
  X_tilde = X_tilde_adj,
  x_split_method = "split_4quantile_x",
  probs = NULL,
  pts = NULL,
  n_group = 30,
  funct_YU = fun_Y_3_int,
  bandwidth = 0.2,
  smooth_curve = TRUE,
  num = 3
)
dev.off()

pdf("splitting_n=5000_setting_A_senario3.pdf",
    width = 9,
    height = 7)
ggplot_posterior_splitting(
  data_sim_3,
  model_adj3,
  X_tilde = X_tilde_adj,
  x_split_method = "split_4quantile_x",
  probs = NULL,
  pts = NULL,
  n_group = 30,
  funct_YU = fun_Y_3_int,
  bandwidth = 0.2,
  smooth_curve = FALSE,
  num = 3
)
dev.off()


n = 6250
data_sim_4 = simulation_data_4(n)
n = dim(data_sim_4)[1]
X_tilde_adj = as.matrix(cbind(rep(1, n), data_sim_4$X, data_sim_4$Z))
X_split_method = "split_4quantile"

model_adj4 <- DDP_ADJ(
  s_seed = 10,
  X_tilde = X_tilde_adj,
  X_W_reg = X_tilde_adj,
  X_split_method,
  probs = NULL,
  pts = NULL,
  R = 4000,
  R_burnin = 3000,
  n_group = 30,
  n = n,
  data_sim_4
)

pdf(
  "splitting_smoothing_n=5000_setting_A_senario4.pdf",
  width = 9,
  height = 7
)
ggplot_posterior_splitting(
  data_sim_4,
  model_adj4,
  X_tilde = X_tilde_adj,
  x_split_method = "split_4quantile_x",
  probs = NULL,
  pts = NULL,
  n_group = 30,
  funct_YU = fun_Y_4_int,
  bandwidth = 0.2,
  smooth_curve = TRUE,
  num = 4
)
dev.off()

pdf("splitting_n=5000_setting_A_senario4.pdf",
    width = 9,
    height = 7)
ggplot_posterior_splitting(
  data_sim_4,
  model_adj4,
  X_tilde = X_tilde_adj,
  x_split_method = "split_4quantile_x",
  probs = NULL,
  pts = NULL,
  n_group = 30,
  funct_YU = fun_Y_4_int,
  bandwidth = 0.2,
  smooth_curve = FALSE,
  num = 4
)
dev.off()
