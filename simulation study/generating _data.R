###################################################################################
# ---              SIMULATION STUDY             ---
# ---                4 settings                 ---
# ---       300 samples x 5000 sample size      ---
###################################################################################

# library 
library(mvtnorm)
library(invgamma)
library(BNPmix)
library(truncnorm)
library(parallel)

###################################################################################
# continuous functions for Y (outcome) 
# given treatment (X) and U (unmeasured confounder)

# --- MIXTURE OF LINEARS ---
fun_Y_1 <- function(X,U){
  return(I(X<5.5)*(+1+2*X+2*U)+I(X>=5.5)*(-16+5*X+2.5*U))
}
fun_Y_1_int <- function(X,U){
  return(fun_Y_1(X,mean(U)))
}

# --- PARABOLA ---
fun_Y_2 <- function(X,U){
  return(-10+2.2*(X-6)^2+4*U)
}
fun_Y_2_int <- function(X,U){
  return(fun_Y_2(X,mean(U)))
}

# --- S FUNCTION ---
fun_Y_3 <- function(X,U){
  return(1.5*sign(X-5)*(abs(X-5)^(1/2))+1.7*U)
}
fun_Y_3_int <- function(X,U){
  return(fun_Y_3(X,mean(U)))
}

# --- HALF HORIZONTAL PARABOLA ---
fun_Y_4 <- function(X,U){
  return(-2*exp(-1.4*(X-6))+0.8*exp(U))
}
fun_Y_4_int <- function(X,U){
  return(-2*exp(-1.4*(X-6))+0.8*mean(exp((U))))
}

###################################################################################
# general informations

n=6250       #number of units for each sample, corrispond to 5000 points 
#(I will cut the tails because too scattered)
sample=300   #number of samples  

###################################################################################

# --- generation of samplers for the 4 settings ---

# generating functions
simulation_data_1<-function(n){
  U=rnorm(n,1,0.2)
  Z=rnorm(n,-1+1.5*U,0.2)
  W=rnorm(n,1-2*U,0.2)
  X=1.5+4*U+rnorm(n,0,0.2)
  Y=rnorm(n, fun_Y_1(X,U),0.3)
  a=cbind(U,W,Z,X,Y)
  a=a[X>quantile(X,0.1) & X<quantile(X,0.9), ]
  return(as.data.frame(a))
}

simulation_data_2<-function(n){
  U=rnorm(n,1,0.2)
  Z=rnorm(n,-1+1.5*U,0.2)
  W=rnorm(n,1-2*U,0.2)
  X=2.5+4*U+rnorm(n,0,0.2)
  Y=rnorm(n, fun_Y_2(X,U),0.2)
  a=cbind(U,W,Z,X,Y)
  a=a[X>quantile(X,0.1) & X<quantile(X,0.9), ]
  return(as.data.frame(a))
}

simulation_data_3<-function(n){
  U=rnorm(n,1,0.2)
  Z=rnorm(n,-1+1.5*U,0.2)
  W=rnorm(n,1-2*U,0.2)
  X=1+4*U+rnorm(n,0,0.2)
  Y=rnorm(n, fun_Y_3(X,U),0.05)
  a=cbind(U,W,Z,X,Y)
  a=a[X>quantile(X,0.1) & X<quantile(X,0.9), ]
  return(as.data.frame(a))
}

simulation_data_4<-function(n){
  U=rnorm(n,1,0.2)
  Z=rnorm(n,-1+1.5*U,0.2)
  W=rnorm(n,1-2*U,0.2)
  X=2.5+4*U+rnorm(n,0,0.2)
  Y=rnorm(n, fun_Y_4(X,U),0.2)
  a=cbind(U,W,Z,X,Y)
  a=a[X>quantile(X,0.1) & X<quantile(X,0.9), ]
  return(as.data.frame(a))
}

# generate data
data_sim_1=lapply(1:sample, function(i) simulation_data_1(n))
data_sim_2=lapply(1:sample, function(i) simulation_data_2(n))
data_sim_3=lapply(1:sample, function(i) simulation_data_3(n))
data_sim_4=lapply(1:sample, function(i) simulation_data_4(n))

save.image("data_simulations.RData")
