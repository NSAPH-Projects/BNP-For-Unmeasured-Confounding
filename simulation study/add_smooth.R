load("results_model_PAPER.RData")

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
  points_y_smooth_02=matrix(NA,nrow=sample,ncol=length(points_x_i))
  points_y_smooth_03=matrix(NA,nrow=sample,ncol=length(points_x_i))
  points_y_smooth_05=matrix(NA,nrow=sample,ncol=length(points_x_i))
  for(i in 1:sample){
    points_y_samples[i,]=sapply(points_x_i, function(x) 
      curve_chains_with(x,data_sim[[i]],post_chian[[i]],split_function_x))
    points_y_samples_NO[i,]=sapply(points_x_i, function(x) 
      curve_chains_without(x,data_sim[[i]],post_chian[[i]],split_function_x))
    points_y_smooth_02[i,]=ksmooth(points_x_i, points_y_samples[i,], "normal", bandwidth = 0.2)$y
    points_y_smooth_03[i,]=ksmooth(points_x_i, points_y_samples[i,], "normal", bandwidth = 0.3)$y
    points_y_smooth_05[i,]=ksmooth(points_x_i, points_y_samples[i,], "normal", bandwidth = 0.5)$y
  }
  
  return(list(x_points=points_x_i,with=points_y_samples,
              without=points_y_samples_NO,
              smooth02=points_y_smooth_02,smooth03=points_y_smooth_03,
              smooth05=points_y_smooth_05))
}

# --- setting 1: ---
start.time <- Sys.time()
median_1s_4q=causal_effect_median(data_sim=data_sim_1, 
                                  post_chian=post_1s_4,
                                  split_function_x=split_4quantile_x)
end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken
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

save.image("results_model_PAPER_smooth.RData")
