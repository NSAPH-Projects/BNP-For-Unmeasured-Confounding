###################################################################################
# ---          SIMULATION FOR THE PAPER         ---
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
# generation of samplers for the 4 settings

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

data_sim_1=lapply(1:sample, function(i) simulation_data_1(n))
data_sim_2=lapply(1:sample, function(i) simulation_data_2(n))
data_sim_3=lapply(1:sample, function(i) simulation_data_3(n))
data_sim_4=lapply(1:sample, function(i) simulation_data_4(n))

###################################################################################
# ---    GIBBS for Dependent Dirichlet Mixture Model      ----
###################################################################################

# basic setting for Gibbs

R=2000            # iteartions
R_burnin=1000     # burn-in
n_group=10        # max number of groups

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

split_6quantile<-function(x_matrix){
  seq_X=c(min(x_matrix),quantile(x_matrix,c(0.1667,0.333,0.5,0.667,0.833)), max(x_matrix))
  X_weights=cbind(x_matrix*I(x_matrix<seq_X[2]),
                  x_matrix*I(x_matrix<seq_X[3] & x_matrix>seq_X[2]),
                  x_matrix*I(x_matrix<seq_X[4] & x_matrix>seq_X[3]),
                  x_matrix*I(x_matrix<seq_X[5] & x_matrix>seq_X[4]),
                  x_matrix*I(x_matrix<seq_X[6] & x_matrix>seq_X[5]),
                  x_matrix*I(x_matrix>seq_X[6]))
  X_weights=(X_weights!=0)*(X_weights)+
    (X_weights==0)*matrix(rnorm(prod(dim(X_weights)),0,0.003),
                          ncol=dim(X_weights)[2],nrow=dim(X_weights)[1])
  
  return(X_weights)
}

###################################################################################

# Gibbs samplers

DDP_Y_4<-function(c, data_sim, split_function){
  
  data_sim=data_sim[[c]]
  
  n=length(data_sim$X)
  
  X_tilde=cbind(rep(1,n),data_sim$X,data_sim$Z)
  X_weights=split_function(data_sim$X)
  
  #prior
  alpha_mu=c(0,0,0)
  alpha_sigma=1
  delta_mu=c(0,0,0)
  delta_sigma=1
  gamma_y_prior=c(1,0.5)
  gamma_w_prior=c(2,2)
  eta_prior=c(0,5)
  
  #initialization
  alpha=matrix(rep(rmvnorm(1,alpha_mu*rep(1,3),alpha_sigma*diag(3)),n_group),
               ncol=n_group, nrow=3)
  delta=rnorm(3,delta_mu,delta_sigma)
  gamma_y=rep(rinvgamma(1,gamma_y_prior[1],gamma_y_prior[2]),n_group)
  gamma_w=rinvgamma(1,gamma_w_prior[1],gamma_w_prior[2])
  eta=matrix(rnorm((n_group-1)*2,eta_prior[1],eta_prior[2]),
             ncol=n_group-1,nrow=dim(X_weights)[2])
  eta_X=pnorm(cbind(X_weights%*%(eta[,1]),X_weights%*%(eta[,2])))
  for (g_k in 3:(n_group-1)) {
    eta_X=cbind(eta_X,pnorm(X_weights%*%(eta[,g_k])))
  }
  eta_X=cbind(eta_X,1)
  tau_prov=cbind(eta_X[,1]*dnorm(data_sim$Y,X_tilde%*%(alpha[,1]),sqrt(gamma_y[1])),
                 eta_X[,2]*(1-eta_X[,1])*dnorm(data_sim$Y,X_tilde%*%(alpha[,2]),sqrt(gamma_y[2])))
  for (g_k in 3:n_group) {
    tau_prov=cbind(tau_prov,eta_X[,g_k]*apply(1-eta_X[,1:(g_k-1)], 1, prod)*
                     dnorm(data_sim$Y,X_tilde%*%(alpha[,g_k]),sqrt(gamma_y[g_k])) )
  }
  tau=tau_prov/apply(tau_prov, 1, sum)
  
  K=sample(1:n_group,n,replace=TRUE)
  table_Ki=rep(0,n_group)
  table_Ki[sort(unique(K))]=table(K)
  
  Q=matrix(NA,nrow=n, ncol=n_group-1)
  beta_x=rep(0,n_group)
  
  #chains
  post_alpha=matrix(NA,nrow=n_group*3, ncol=R)
  post_delta=matrix(NA,nrow=3, ncol=R)
  post_gamma_y=matrix(NA,nrow=n_group, ncol=R)
  post_gamma_w=rep(NA,R)
  K_all=matrix(NA,nrow=n, ncol=R)
  post_beta_x=matrix(NA,nrow=n_group, ncol=R)
  post_intercept=matrix(NA,nrow=n_group, ncol=R)
  post_eta=matrix(NA,nrow=(n_group-1)*(dim(X_weights)[2]), ncol=R)
  
  for (r in 1:R){
    
    #delta  ( parameters of the mean in the distribution of W|X,Z )
    V_delta=diag(3)/delta_sigma+t(X_tilde)%*%X_tilde/gamma_w
    m_delta=delta_mu/delta_sigma+t(X_tilde)%*%data_sim$W/gamma_w
    delta=rmvnorm(1,solve(V_delta)%*%m_delta, solve(V_delta))
    
    #gamma_w    ( parameter of the variance in the distribution of W|X,Z )
    g_w=sum((data_sim$W-X_tilde%*%t(delta))^2)
    gamma_w=rinvgamma(1,gamma_w_prior[1]+n/2,gamma_w_prior[2]+g_w/2)
    
    #for compute eta (regression parameters in the weights)
    # we need to compute mu_Q and Q
    
    #mu_Q
    mu_Q=cbind((tau[,1]),(tau[,2]/(1-tau[,1])))
    if (n_group>3){
      mu_Q=cbind(mu_Q,sapply(3:(n_group-1), function(l) (tau[,l]/(1-apply(tau[,1:(l-1)],1,sum)))))
    }
    # for resolve problems in the bounds
    mu_Q[which(is.nan(mu_Q))]=1
    mu_Q[which(mu_Q>1)]=1
    mu_Q=mu_Q-9.9e-15*(mu_Q>(1-1e-16))
    mu_Q=mu_Q+1e-16*(mu_Q<1e-16)
    
    #Q
    for (i in 1:n){
      for (l in 1:(min(K[i],n_group-1))) {
        if (l<K[i]){
          Q[i,l]=rtruncnorm(1,b=0,mean=qnorm(mu_Q[i,l]))
        }else{
          Q[i,l]=rtruncnorm(1,a=0,mean=qnorm(mu_Q[i,l]))
        }
      }
    }
    
    #eta
    check_q=sapply(1:(n_group-1), function(l) length(Q[K>=l,l]))
    v_eta=lapply(which(check_q>2), function(l) 1/eta_prior[2]+t(X_weights[K>=l,])%*%X_weights[K>=l,])
    m_eta=sapply(which(check_q>2), function(l) eta_prior[1]/eta_prior[2]+
                   t(X_weights[K>=l,])%*%Q[K>=l,l]/gamma_y[l])
    
    if (any(check_q<3)){
      for(l in which(check_q==1 | check_q==2)){
        v_eta[[l]]=diag(dim(X_weights)[2])/eta_prior[2]
        m_eta=cbind(m_eta,rep(eta_prior[1]/eta_prior[2],2))
      }
    }
    eta[,which(table_Ki[-n_group]>0)]=sapply(which(table_Ki[-n_group]>0), 
                                             function(l) rmvnorm(1,solve(v_eta[[l]])%*%m_eta[,l]))
    eta[,which(table_Ki[-n_group]==0)]=rnorm(length(which(table_Ki[-n_group]==0))*2,eta_prior[1],sqrt(eta_prior[2]))
    
    
    #tau ( recursive weigths in the mixture )
    eta_X=pnorm(cbind(X_weights%*%(eta[,1]),X_weights%*%(eta[,2])))
    for (g_k in 3:(n_group-1)) {
      eta_X=cbind(eta_X,pnorm(X_weights%*%(eta[,g_k])))
    }
    eta_X=cbind(eta_X,1)
    
    tau_prov=cbind(eta_X[,1]*dnorm(data_sim$Y,X_tilde%*%(alpha[,1]),sqrt(gamma_y[1])),
                   eta_X[,2]*(1-eta_X[,1])*dnorm(data_sim$Y,X_tilde%*%(alpha[,2]),sqrt(gamma_y[2])))
    for (g_k in 3:n_group) {
      tau_prov=cbind(tau_prov,eta_X[,g_k]*apply(1-eta_X[,1:(g_k-1)], 1, prod)*
                       dnorm(data_sim$Y,X_tilde%*%(alpha[,g_k]),sqrt(gamma_y[g_k])) )
    }
    
    tau=tau_prov/apply(tau_prov, 1, sum)
    for (i in which(is.na(tau[,1]))){
      tau[i,]=1/n_group
    }
    
    #K
    K=sapply(1:n, function(i) (1:n_group)%*%rmultinom(1,1,tau[i,]))
    table_Ki=rep(0,n_group)
    table_Ki[sort(unique(K))]=table(K)
    
    #alpha + gamma_y ( parameters of the mean + variance in each cluster of the distribution of Y|X,Z )
    for (g_k in 1:n_group){
      if(length(which(K==g_k))==1){
        V=diag(3)/alpha_sigma+t(t(X_tilde[K==g_k,]))%*%t(X_tilde[K==g_k,])/gamma_y[g_k]
        M=alpha_mu/alpha_sigma+t(t(X_tilde[K==g_k,]))%*%t(data_sim$Y[K==g_k])/gamma_y[g_k]
      }else{
        V=diag(3)/alpha_sigma+t(X_tilde[K==g_k,])%*%X_tilde[K==g_k,]/gamma_y[g_k]
        M=alpha_mu/alpha_sigma+t(X_tilde[K==g_k,])%*%data_sim$Y[K==g_k]/gamma_y[g_k]
      }
      alpha[,g_k]=t(rmvnorm(1,solve(V)%*%M, solve(V)))
      
      G=sum((data_sim$Y[K==g_k]-X_tilde[K==g_k,]%*%(alpha[,g_k]))^2)
      gamma_y[g_k]=rinvgamma(1,gamma_y_prior[1]+(sum(K==g_k))/2,gamma_y_prior[2]+G/2)
    }
    
    
    #beta_x
    beta_x=alpha[2,]-alpha[3,]*delta[1,2]/delta[1,3]
    intercept=sapply(1:n_group, function(g) 
      mean(alpha[1,g]+alpha[3,g]*mean(data_sim$Z)+alpha[3,g]*delta[1,2]/delta[1,3]*mean(data_sim$X)))
    
    # ---------- save the chains  ------
    post_alpha[,r]=c(alpha)
    post_delta[,r]=delta[,1:3]
    post_gamma_y[,r]=gamma_y
    post_gamma_w[r]=gamma_w
    K_all[,r]=K
    post_beta_x[,r]=beta_x
    post_intercept[,r]=intercept
    post_eta[,r]=c(eta)
    
    #check
    if (r%%500==0) print(r)
  }
  
  print(c)
  
  return(list(post_alpha=post_alpha[,(R_burnin+1):R],
              post_beta_x=post_beta_x[,(R_burnin+1):R], 
              post_intercept=post_intercept[,(R_burnin+1):R],
              post_eta=post_eta[,(R_burnin+1):R]))
}

DDP_Y_6<-function(c, data_sim, split_function){
  
  data_sim=data_sim[[c]]
  
  n=length(data_sim$X)
  
  X_tilde=cbind(rep(1,n),data_sim$X,data_sim$Z)
  X_weights=split_function(data_sim$X)
  
  #prior
  alpha_mu=c(0,0,0)
  alpha_sigma=1
  delta_mu=c(0,0,0)
  delta_sigma=1
  gamma_y_prior=c(1,0.5)
  gamma_w_prior=c(2,2)
  eta_prior=c(0,5)
  
  #initialization
  alpha=matrix(rep(rmvnorm(1,alpha_mu*rep(1,3),alpha_sigma*diag(3)),n_group),
               ncol=n_group, nrow=3)
  delta=rnorm(3,delta_mu,delta_sigma)
  gamma_y=rep(rinvgamma(1,gamma_y_prior[1],gamma_y_prior[2]),n_group)
  gamma_w=rinvgamma(1,gamma_w_prior[1],gamma_w_prior[2])
  eta=matrix(rnorm((n_group-1)*2,eta_prior[1],eta_prior[2]),
             ncol=n_group-1,nrow=dim(X_weights)[2])
  eta_X=pnorm(cbind(X_weights%*%(eta[,1]),X_weights%*%(eta[,2])))
  for (g_k in 3:(n_group-1)) {
    eta_X=cbind(eta_X,pnorm(X_weights%*%(eta[,g_k])))
  }
  eta_X=cbind(eta_X,1)
  tau_prov=cbind(eta_X[,1]*dnorm(data_sim$Y,X_tilde%*%(alpha[,1]),sqrt(gamma_y[1])),
                 eta_X[,2]*(1-eta_X[,1])*dnorm(data_sim$Y,X_tilde%*%(alpha[,2]),sqrt(gamma_y[2])))
  for (g_k in 3:n_group) {
    tau_prov=cbind(tau_prov,eta_X[,g_k]*apply(1-eta_X[,1:(g_k-1)], 1, prod)*
                     dnorm(data_sim$Y,X_tilde%*%(alpha[,g_k]),sqrt(gamma_y[g_k])) )
  }
  tau=tau_prov/apply(tau_prov, 1, sum)
  
  K=sample(1:n_group,n,replace=TRUE)
  table_Ki=rep(0,n_group)
  table_Ki[sort(unique(K))]=table(K)
  
  Q=matrix(NA,nrow=n, ncol=n_group-1)
  beta_x=rep(0,n_group)
  
  #chains
  post_alpha=matrix(NA,nrow=n_group*3, ncol=R)
  post_delta=matrix(NA,nrow=3, ncol=R)
  post_gamma_y=matrix(NA,nrow=n_group, ncol=R)
  post_gamma_w=rep(NA,R)
  K_all=matrix(NA,nrow=n, ncol=R)
  post_beta_x=matrix(NA,nrow=n_group, ncol=R)
  post_intercept=matrix(NA,nrow=n_group, ncol=R)
  post_eta=matrix(NA,nrow=(n_group-1)*(dim(X_weights)[2]), ncol=R)
  
  for (r in 1:R){
    
    #delta  ( parameters of the mean in the distribution of W|X,Z )
    V_delta=diag(3)/delta_sigma+t(X_tilde)%*%X_tilde/gamma_w
    m_delta=delta_mu/delta_sigma+t(X_tilde)%*%data_sim$W/gamma_w
    delta=rmvnorm(1,solve(V_delta)%*%m_delta, solve(V_delta))
    
    #gamma_w    ( parameter of the variance in the distribution of W|X,Z )
    g_w=sum((data_sim$W-X_tilde%*%t(delta))^2)
    gamma_w=rinvgamma(1,gamma_w_prior[1]+n/2,gamma_w_prior[2]+g_w/2)
    
    #for compute eta (regression parameters in the weights)
    # we need to compute mu_Q and Q
    
    #mu_Q
    mu_Q=cbind((tau[,1]),(tau[,2]/(1-tau[,1])))
    if (n_group>3){
      mu_Q=cbind(mu_Q,sapply(3:(n_group-1), function(l) (tau[,l]/(1-apply(tau[,1:(l-1)],1,sum)))))
    }
    # for resolve problems in the bounds
    mu_Q[which(is.nan(mu_Q))]=1
    mu_Q[which(mu_Q>1)]=1
    mu_Q=mu_Q-9.9e-15*(mu_Q>(1-1e-16))
    mu_Q=mu_Q+1e-16*(mu_Q<1e-16)
    
    #Q
    for (i in 1:n){
      for (l in 1:(min(K[i],n_group-1))) {
        if (l<K[i]){
          Q[i,l]=rtruncnorm(1,b=0,mean=qnorm(mu_Q[i,l]))
        }else{
          Q[i,l]=rtruncnorm(1,a=0,mean=qnorm(mu_Q[i,l]))
        }
      }
    }
    
    #eta
    check_q=sapply(1:(n_group-1), function(l) length(Q[K>=l,l]))
    v_eta=lapply(which(check_q>4), function(l) 1/eta_prior[2]+t(X_weights[K>=l,])%*%X_weights[K>=l,])
    m_eta=sapply(which(check_q>4), function(l) eta_prior[1]/eta_prior[2]+
                   t(X_weights[K>=l,])%*%Q[K>=l,l]/gamma_y[l])
    
    if (any(check_q<5)){
      for(l in which(check_q==1 | check_q==2 | check_q==3 | check_q==4)){
        v_eta[[l]]=diag(dim(X_weights)[2])/eta_prior[2]
        m_eta=cbind(m_eta,rep(eta_prior[1]/eta_prior[2],2))
      }
    }
    eta[,which(table_Ki[-n_group]>0)]=sapply(which(table_Ki[-n_group]>0), 
                                             function(l) rmvnorm(1,solve(v_eta[[l]])%*%m_eta[,l]))
    eta[,which(table_Ki[-n_group]==0)]=rnorm(length(which(table_Ki[-n_group]==0))*2,eta_prior[1],sqrt(eta_prior[2]))
    
    
    #tau ( recursive weigths in the mixture )
    eta_X=pnorm(cbind(X_weights%*%(eta[,1]),X_weights%*%(eta[,2])))
    for (g_k in 3:(n_group-1)) {
      eta_X=cbind(eta_X,pnorm(X_weights%*%(eta[,g_k])))
    }
    eta_X=cbind(eta_X,1)
    
    tau_prov=cbind(eta_X[,1]*dnorm(data_sim$Y,X_tilde%*%(alpha[,1]),sqrt(gamma_y[1])),
                   eta_X[,2]*(1-eta_X[,1])*dnorm(data_sim$Y,X_tilde%*%(alpha[,2]),sqrt(gamma_y[2])))
    for (g_k in 3:n_group) {
      tau_prov=cbind(tau_prov,eta_X[,g_k]*apply(1-eta_X[,1:(g_k-1)], 1, prod)*
                       dnorm(data_sim$Y,X_tilde%*%(alpha[,g_k]),sqrt(gamma_y[g_k])) )
    }
    
    tau=tau_prov/apply(tau_prov, 1, sum)
    for (i in which(is.na(tau[,1]))){
      tau[i,]=1/n_group
    }
    
    #K
    K=sapply(1:n, function(i) (1:n_group)%*%rmultinom(1,1,tau[i,]))
    table_Ki=rep(0,n_group)
    table_Ki[sort(unique(K))]=table(K)
    
    #alpha + gamma_y ( parameters of the mean + variance in each cluster of the distribution of Y|X,Z )
    for (g_k in 1:n_group){
      if(length(which(K==g_k))==1){
        V=diag(3)/alpha_sigma+t(t(X_tilde[K==g_k,]))%*%t(X_tilde[K==g_k,])/gamma_y[g_k]
        M=alpha_mu/alpha_sigma+t(t(X_tilde[K==g_k,]))%*%t(data_sim$Y[K==g_k])/gamma_y[g_k]
      }else{
        V=diag(3)/alpha_sigma+t(X_tilde[K==g_k,])%*%X_tilde[K==g_k,]/gamma_y[g_k]
        M=alpha_mu/alpha_sigma+t(X_tilde[K==g_k,])%*%data_sim$Y[K==g_k]/gamma_y[g_k]
      }
      alpha[,g_k]=t(rmvnorm(1,solve(V)%*%M, solve(V)))
      
      G=sum((data_sim$Y[K==g_k]-X_tilde[K==g_k,]%*%(alpha[,g_k]))^2)
      gamma_y[g_k]=rinvgamma(1,gamma_y_prior[1]+(sum(K==g_k))/2,gamma_y_prior[2]+G/2)
    }
    
    
    #beta_x
    beta_x=alpha[2,]-alpha[3,]*delta[1,2]/delta[1,3]
    intercept=sapply(1:n_group, function(g) 
      mean(alpha[1,g]+alpha[3,g]*mean(data_sim$Z)+alpha[3,g]*delta[1,2]/delta[1,3]*mean(data_sim$X)))
    
    # ---------- save the chains  ------
    post_alpha[,r]=c(alpha)
    post_delta[,r]=delta[,1:3]
    post_gamma_y[,r]=gamma_y
    post_gamma_w[r]=gamma_w
    K_all[,r]=K
    post_beta_x[,r]=beta_x
    post_intercept[,r]=intercept
    post_eta[,r]=c(eta)
    
    #check
    if (r%%500==0) print(r)
  }
  
  print(c)
  
  return(list(post_alpha=post_alpha[,(R_burnin+1):R],
              post_beta_x=post_beta_x[,(R_burnin+1):R], 
              post_intercept=post_intercept[,(R_burnin+1):R],
              post_eta=post_eta[,(R_burnin+1):R]))
}

###################################################################################

# --- setting 1: ---
start.time <- Sys.time()
post_1s_4=mclapply(1:sample, DDP_Y_4, data_sim=data_sim_1, 
                         split_function=split_4quantile, mc.cores=6)
end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken

start.time <- Sys.time()
post_1s_6=mclapply(1:sample, DDP_Y_6, data_sim=data_sim_1, 
                         split_function=split_6quantile, mc.cores=6)
end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken

# --- setting 2: ---
start.time <- Sys.time()
post_2s_4=mclapply(1:sample, DDP_Y_4, data_sim=data_sim_2, 
                         split_function=split_4quantile, mc.cores=6)
end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken

start.time <- Sys.time()
post_2s_6=mclapply(1:sample, DDP_Y_6, data_sim=data_sim_2, 
                         split_function=split_6quantile, mc.cores=6)
end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken

# --- setting 3: ---
start.time <- Sys.time()
post_3s_4=mclapply(1:sample, DDP_Y_4, data_sim=data_sim_3, 
                         split_function=split_4quantile, mc.cores=6)
end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken

start.time <- Sys.time()
post_3s_6=mclapply(1:sample, DDP_Y_6, data_sim=data_sim_3, 
                         split_function=split_6quantile, mc.cores=6)
end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken

# --- setting 4: ---
start.time <- Sys.time()
post_4s_4=mclapply(1:sample, DDP_Y_4, data_sim=data_sim_4, 
                         split_function=split_4quantile, mc.cores=6)
end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken

start.time <- Sys.time()
post_4s_6=mclapply(1:sample, DDP_Y_6, data_sim=data_sim_4, 
                         split_function=split_6quantile, mc.cores=6)
end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken

save.image("results_model_PAPER.RData")

###################################################################################
# ---    compute MEDIANS for the posterior distributions      ----
###################################################################################


split_4quantile_x<-function(x_matrix,x){
  seq_X=c(min(x_matrix),quantile(x_matrix,c(0.25,0.5,0.75)), max(x_matrix))
  X_weights=cbind(x*I(x<seq_X[2]),
                  x*I(x<seq_X[3] & x>seq_X[2]),
                  x*I(x<seq_X[4] & x>seq_X[3]),
                  x*I(x>seq_X[4]))
  X_weights=(X_weights!=0)*(X_weights)+
    (X_weights==0)*matrix(rnorm(prod(dim(X_weights)),0,0.003),
                          ncol=dim(X_weights)[2],nrow=dim(X_weights)[1])
  
  return(X_weights)
}

split_6quantile_x<-function(x_matrix,x){
  seq_X=c(min(x_matrix),quantile(x_matrix,c(0.1667,0.333,0.5,0.667,0.833)), max(x_matrix))
  X_weights=cbind(x*I(x<seq_X[2]),
                  x*I(x<seq_X[3] & x>seq_X[2]),
                  x*I(x<seq_X[4] & x>seq_X[3]),
                  x*I(x<seq_X[5] & x>seq_X[4]),
                  x*I(x<seq_X[6] & x>seq_X[5]),
                  x*I(x>seq_X[6]))
  X_weights=(X_weights!=0)*(X_weights)+
    (X_weights==0)*matrix(rnorm(prod(dim(X_weights)),0,0.003),
                          ncol=dim(X_weights)[2],nrow=dim(X_weights)[1])
  
  return(X_weights)
}

#################################################################################

curve_chains_with<- function(x,data_sim,post_chian,split_function_x){
  
  X_weights=split_function_x(data_sim$X,x)
  
  c_eta_X=pnorm(cbind(t(X_weights%*%(post_chian$post_eta[1:(dim(X_weights)[2]),])),
                      t(X_weights%*%(post_chian$post_eta[(dim(X_weights)[2]+1):(dim(X_weights)[2]*2),]))))
  for (g in 3:(n_group-1)) {
    c_eta_X=cbind(c_eta_X,t(pnorm(X_weights%*%(post_chian$post_eta[((g-1)*(dim(X_weights)[2])+1):(g*(dim(X_weights)[2])),]))))
  }
  c_eta_X=cbind(c_eta_X,1)
  
  t=cbind(c_eta_X[,1], c_eta_X[,2]*(1-c_eta_X[,1]),
          sapply(3:n_group, function(g) c_eta_X[,g]*apply(1-c_eta_X[,1:(g-1)], 1, prod)))
  
  clusters=apply(t, 1, which.max)
  values=sapply(1:length(clusters), function(c) post_chian$post_intercept[clusters[c],c]+
                  x%*%post_chian$post_beta_x[clusters[c],c])
  
  return(quantile(values, prob=0.5, na.rm=TRUE))
}

curve_chains_without<- function(x,data_sim,post_chian,split_function_x){
  
  X_weights=split_function_x(data_sim$X,x)
  
  c_eta_X=pnorm(cbind(t(X_weights%*%(post_chian$post_eta[1:(dim(X_weights)[2]),])),
                      t(X_weights%*%(post_chian$post_eta[(dim(X_weights)[2]+1):(dim(X_weights)[2]*2),]))))
  for (g in 3:(n_group-1)) {
    c_eta_X=cbind(c_eta_X,t(pnorm(X_weights%*%(post_chian$post_eta[((g-1)*(dim(X_weights)[2])+1):(g*(dim(X_weights)[2])),]))))
  }
  c_eta_X=cbind(c_eta_X,1)
  
  t=cbind(c_eta_X[,1], c_eta_X[,2]*(1-c_eta_X[,1]),
          sapply(3:n_group, function(g) c_eta_X[,g]*apply(1-c_eta_X[,1:(g-1)], 1, prod)))
  
  clusters=apply(t, 1, which.max)
  values=sapply(1:length(clusters), function(c) 
    post_chian$post_alpha[clusters[c]*3-2,c]+
      post_chian$post_alpha[clusters[c]*3,c]*mean(data_sim$Z)+
      x%*%post_chian$post_alpha[clusters[c]*3-1,c])
  
  return(quantile(values, prob=0.5, na.rm=TRUE))
}

causal_effect_median<- function(data_sim,post_chian,split_function_x){
  
  min_x=mean(unlist(sapply(1:sample, function(i) min(data_sim[[i]]$X))))
  max_x=mean(unlist(sapply(1:sample, function(i) max(data_sim[[i]]$X))))
  
  # causal effect for a grid of points
  points_x_i=seq(min_x,max_x,length.out=200)
  points_y_samples=matrix(NA,nrow=sample,ncol=length(points_x_i))
  points_y_samples_NO=matrix(NA,nrow=sample,ncol=length(points_x_i))
  points_y_smooth=matrix(NA,nrow=sample,ncol=length(points_x_i))
  for(i in 1:sample){
    points_y_samples[i,]=sapply(points_x_i, function(x) 
      curve_chains_with(x,data_sim[[i]],post_chian[[i]],split_function_x))
    points_y_samples_NO[i,]=sapply(points_x_i, function(x) 
      curve_chains_without(x,data_sim[[i]],post_chian[[i]],split_function_x))
    points_y_smooth[i,]=ksmooth(points_x_i, points_y_samples[i,], "normal", bandwidth = 0.1)$y
  }
  
  return(list(x_points=points_x_i,with=points_y_samples,
              without=points_y_samples_NO,smooth=points_y_smooth))
}

# --- setting 1: ---
median_1s_4q=causal_effect_median(data_sim=data_sim_1, 
                                       post_chian=post_1s_4,
                                       split_function_x=split_4quantile_x)
median_1s_6q=causal_effect_median(data_sim=data_sim_1, 
                                       post_chian=post_1s_6,
                                       split_function_x=split_6quantile_x)

# --- setting 2: ---
median_2s_4q=causal_effect_median(data_sim=data_sim_2, 
                                       post_chian=post_2s_4,
                                       split_function_x=split_4quantile_x)
median_2s_6q=causal_effect_median(data_sim=data_sim_2, 
                                       post_chian=post_2s_6,
                                       split_function_x=split_6quantile_x)

# --- setting 3: ---
median_3s_4q=causal_effect_median(data_sim=data_sim_3, 
                                       post_chian=post_3s_4,
                                       split_function_x=split_4quantile_x)
median_3s_6q=causal_effect_median(data_sim=data_sim_3, 
                                       post_chian=post_3s_6,
                                       split_function_x=split_6quantile_x)

# --- setting 4: ---
median_4s_4q=causal_effect_median(data_sim=data_sim_4, 
                                       post_chian=post_4s_4,
                                       split_function_x=split_4quantile_x)
median_4s_6q=causal_effect_median(data_sim=data_sim_4, 
                                       post_chian=post_4s_6,
                                       split_function_x=split_6quantile_x)

save.image("results_model_PAPER.RData")
