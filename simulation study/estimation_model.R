###################################################################################
# ---              SIMULATION STUDY             ---
# ---                4 settings                 ---
###################################################################################

# library 
library(mvtnorm)
library(invgamma)
library(BNPmix)
library(truncnorm)
library(parallel)

# upload data
load("data_simulations.RData")

# upload Gibbs sampler
source()

###################################################################################

# basic setting for Gibbs

R=2000            # iteartions
R_burnin=1000     # burn-in
n_group=10        # max number of groups

###################################################################################
# --- SECTION: SIMULATION STUDY --- 

# --- setting 1: ---
post_1s_4=mclapply(1:sample, DDP_Y_4, data_sim=data_sim_1, 
                   split_function=split_4quantile, mc.cores=6)
# --- setting 2: ---
post_2s_4=mclapply(1:sample, DDP_Y_4, data_sim=data_sim_2, 
                   split_function=split_4quantile, mc.cores=6)
# --- setting 3: ---
post_3s_4=mclapply(1:sample, DDP_Y_4, data_sim=data_sim_3, 
                   split_function=split_4quantile, mc.cores=6)
# --- setting 4: ---
post_4s_4=mclapply(1:sample, DDP_Y_4, data_sim=data_sim_4, 
                   split_function=split_4quantile, mc.cores=6)

###################################################################################
# --- APPENDIX ----

# --- setting 1: ---
post_1s_6=mclapply(1:sample, DDP_Y_6, data_sim=data_sim_1, 
                   split_function=split_6quantile, mc.cores=6)

# --- setting 2: ---
post_2s_6=mclapply(1:sample, DDP_Y_6, data_sim=data_sim_2, 
                   split_function=split_6quantile, mc.cores=6)

# --- setting 3: ---
post_3s_6=mclapply(1:sample, DDP_Y_6, data_sim=data_sim_3, 
                   split_function=split_6quantile, mc.cores=6)

# --- setting 4: ---
post_4s_6=mclapply(1:sample, DDP_Y_6, data_sim=data_sim_4, 
                   split_function=split_6quantile, mc.cores=6)


###################################################################################
# ---    compute MEDIANS for the posterior distributions of CERF     ----
###################################################################################

# --- SECTION: SIMULATION STUDY --- 

# --- setting 1: ---
median_1s_4q=causal_effect_median(data_sim=data_sim_1, 
                                  post_chian=post_1s_4,
                                  split_function_x=split_4quantile_x)
# --- setting 2: ---
median_2s_4q=causal_effect_median(data_sim=data_sim_2, 
                                  post_chian=post_2s_4,
                                  split_function_x=split_4quantile_x)
# --- setting 3: ---
median_3s_4q=causal_effect_median(data_sim=data_sim_3, 
                                  post_chian=post_3s_4,
                                  split_function_x=split_4quantile_x)
# --- setting 4: ---
median_4s_4q=causal_effect_median(data_sim=data_sim_4, 
                                  post_chian=post_4s_4,
                                  split_function_x=split_4quantile_x)

###################################################################################
# --- APPENDIX ----

# --- setting 1: ---
median_1s_6q=causal_effect_median(data_sim=data_sim_1, 
                                  post_chian=post_1s_6,
                                  split_function_x=split_6quantile_x)
# --- setting 2: ---
median_2s_6q=causal_effect_median(data_sim=data_sim_2, 
                                  post_chian=post_2s_6,
                                  split_function_x=split_6quantile_x)
# --- setting 3: ---
median_3s_6q=causal_effect_median(data_sim=data_sim_3, 
                                  post_chian=post_3s_6,
                                  split_function_x=split_6quantile_x)
# --- setting 4: ---
median_4s_6q=causal_effect_median(data_sim=data_sim_4, 
                                  post_chian=post_4s_6,
                                  split_function_x=split_6quantile_x)


###################################################################################
# save results 
save.image("results_model.RData")
