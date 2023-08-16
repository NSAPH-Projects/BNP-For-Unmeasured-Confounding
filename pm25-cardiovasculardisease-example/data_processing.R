################################################################################
# ----        analysis dataset       ---
################################################################################

# libraries
rm(list = ls())
library(data.table)

# load cleaning data functions
source('Georgias_code/ReformData_function.R')
source('Georgias_code/DataClean_function.R')

# load datset
load("full_2013_11_12_6M_closetF.RData")

################################################################################

# Specify treatment variable 
trt_col <- 'avgPM'

data_clean <- DataClean(copy(medic_full_dta), trt = trt_col, max_exposure = NULL,
                    drop_missing = TRUE)
data_clean=as.data.frame(data_clean)

################################################################################

# ----  PRELIMINARY ANALYSIS ----

# plots for exposure and outcomes

hist(data_clean$Exposure, main=" ", xlab="PM2.5",nclass=80, freq=FALSE)
legend("topright", legend=paste0(names(quant_exp)," = ", round(quant_exp,2)), 
       title="quantiles")

hist(data_clean$lograte, main=" ", xlab="log cardiovascular diseases rate",
     nclass=80, freq=FALSE)
legend("topright", legend=paste0(names(quant_out)," = ", round(quant_out,2)), 
       title="quantiles")

################################################################################

# ----  CHOICE OF CONFOUNDER ----

conf_var=c(4:19,21:23,25:28)
names_conf=colnames(data_clean)[conf_var]

# correlation with treatment 
cor_conf=sapply(conf_var, function(i) cor(data_clean[,i],data_clean$Exposure))

plot(cor_conf[order(cor_conf)],pch=19, 
     main="correlation  with PM2.5", ylab="cor", xlab=" ", xaxt="n")
axis(1, at=1:(length(names_conf)), labels=names_conf[order(cor_conf)], cex.axis=0.5, las=2)


# correlation with outcome
cor_conf=sapply(conf_var, function(i) cor(data_clean[,i],data_clean$lograte))

par(mfrow=c(1,1))
plot(cor_conf[order(cor_conf)],pch=19, 
     main="correlation  with log rate cvd", ylab="cor", xlab=" ", xaxt="n")
axis(1, at=1:(length(names_conf)), labels=names_conf[order(cor_conf)], cex.axis=0.5, las=2)

# choosing the income
#  logMedHInc00_13

cor(data_clean$logMedHInc00_13,data_clean$Exposure)
cor(data_clean$logMedHInc00_13,data_clean$lograte)

################################################################################

# save cleaned dataset
save(data_clean, file="data_allariables_cleaned.RData")
