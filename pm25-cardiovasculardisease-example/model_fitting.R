################################################################################
# ----        MODEL FITTING       ---
################################################################################

#libraries
library(mvtnorm)
library(invgamma)
library(truncnorm)

#functions
source("src/BNP_NC.R")
source("src/BNP_NoAdj.R")
source("src/CERF.R")
source("src/helper.R")

#data
load("data_allariables_cleaned.RData")

################################################################################

# --- prepare dataset for models functions ---

data_clean$logPctblPvt2013 <- log(data_clean$PctblPvt2013)
data_clean$logPctblPvt2013[data_clean$logPctblPvt2013 == "-Inf"] = -4 

data_analysis = as.data.frame(
  cbind(
    Y = data_clean$lograte,
    X = data_clean$Exposure,
    U = data_clean$logMedHInc00_13,
    Z = data_clean$PctOccupied,
    W = data_clean$PctOwnHs2013
  )
)

# --- fitting models for CERF ---

estimated_CERF=cerf(data_analysis, 
                X_split_method = split_4quantile, 
                x_split_method = split_4quantile_x, 
                n_group = 10, R = 2000, R_burnin = 1000, n=5362)

################################################################################

# save results
save(estimated_CERF, file="echains_APPLICATION.RData")
