# Look at EDA.R for a look in the exposure distribution and the
# plots of standardized difference at different exposure levels.

rm(list = ls())
library(MASS)
library(data.table)
library(mvnfast)
library(truncnorm)
setwd('~/Documents/Research/Completed/Causal_ER_curve/Zipcode_analysis/')
source('functions/ReformData_function.R')
source('functions/DataClean_function.R')

file_source <- list.files('~/Github/LERCA/R/', pattern = '*.R$')
sapply(paste0('~/Github/LERCA/R/', file_source), source, .GlobalEnv)

process <- 1

# ---- Specify -----

trt_col <- 'avgPM'

# ------------------

load('Data/full_2013_11_12_6M_closeF.Rdata')
subdta <- DataClean(copy(medic_full_dta), trt = trt_col, max_exposure = NULL,
                    drop_missing = TRUE)
dta <- ReformData(copy(subdta))

coord_col <- which(names(dta) %in% c('Latitude.zip', 'Longitude.zip'))
er_dta <- dta[, - coord_col, with = FALSE]
er_dta[, zipcode_R := NULL]
out_col <- which(names(er_dta) == 'lograte')
trt_col <- which(names(er_dta) == 'Exposure')
setnames(er_dta, c('Exposure', 'lograte'), c('X', 'Y'))

covs_col <- setdiff(1:ncol(er_dta), c(out_col, trt_col))
conf_names <- names(er_dta)[covs_col]
setnames(er_dta, conf_names, paste0('C', 1 : length(covs_col)))

er_dta <- er_dta[, c(out_col, trt_col, covs_col), with = FALSE]



post_means <- ER$y[, 1, ]  # There is only 1 chain.

chain_waic <- WAIC(lerca = lerca, dta = er_dta)


