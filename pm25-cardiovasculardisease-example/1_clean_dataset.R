################################################################################
# ----        analysis dataset       ---
################################################################################

# libraries
rm(list = ls())
library(data.table)
library(CondIndTests)
setwd('/n/dominici_nsaph_l3/projects/pm25-cardiovasculardisease-bnp/')

# load Giorgia's functions
source('Georgias_code/ReformData_function.R')
source('Georgias_code/DataClean_function.R')

# Specify treatment variable 
trt_col <- 'avgPM'

# load datset
load("full_2013_11_12_6M_closetF.RData")

data_clean <- DataClean(copy(medic_full_dta), trt = trt_col, max_exposure = NULL,
                    drop_missing = TRUE)
data_clean=as.data.frame(data_clean)

# we can clean more the dataset
# scaling and dicotomizing some variables
# I don't think it is necessary for our case
# anyway this is the code:
#dataset_cleaned <- ReformData(copy(dataset_cleaned))

################################################################################

# ----  PRELIMINARY ANALYSIS ----

quant_exp=quantile(data_clean$Exposure)
quant_out=quantile(data_clean$lograte)
cor_exp_out=cor(data_clean$Exposure,data_clean$lograte)

# plots for exposure and outcomes

hist(data_clean$Exposure, main=" ", xlab="PM2.5",nclass=80, freq=FALSE)
legend("topright", legend=paste0(names(quant_exp)," = ", round(quant_exp,2)), 
       title="quantiles")

hist(data_clean$lograte, main=" ", xlab="log cardiovascular diseases rate",
     nclass=80, freq=FALSE)
legend("topright", legend=paste0(names(quant_out)," = ", round(quant_out,2)), 
       title="quantiles")

plot(data_clean$Exposure,data_clean$lograte, 
     pch=19, cex=0.5, xlab="PM2.5", ylab="log CVD rate")
curve_values=ksmooth(data_clean$Exposure,data_clean$lograte, "normal", 
                     bandwidth = 1)
points(curve_values$x,curve_values$y, col="blue",cex=0.3)
legend("topright", legend=paste0("cor = ",round(cor_exp_out,3)))


################################################################################

# ----  CHOICE OF CONFOUNDER ----

conf_var=c(4:19,21:23,25:28)
names_conf=colnames(data_clean)[conf_var]

# correlation with treatment 
cor_conf=sapply(conf_var, function(i) cor(data_clean[,i],data_clean$Exposure))

plot(cor_conf[order(cor_conf)],pch=19, 
     main="correlation  with PM2.5", ylab="cor", xlab=" ", xaxt="n")
axis(1, at=1:(length(names_conf)), labels=names_conf[order(cor_conf)], cex.axis=0.5, las=2)

# zoom in tails
par(mfrow=c(1,2))
plot(cor_conf[order(cor_conf)][1:5], pch=19, 
     main="correlation  with PM2.5 - zoom in left tail", ylab="cor", xlab=" ", xaxt="n")
axis(1, at=1:5, labels=names_conf[order(cor_conf)][1:5], cex.axis=0.6, las=2)

plot(cor_conf[order(cor_conf)][19:23],pch=19, 
     main="correlation  with PM2.5 - zoom in right tail", ylab="cor", xlab=" ", xaxt="n")
axis(1, at=1:5, labels=names_conf[order(cor_conf)][19:23], cex.axis=0.6, las=2)


# correlation with outcome
cor_conf=sapply(conf_var, function(i) cor(data_clean[,i],data_clean$lograte))

par(mfrow=c(1,1))
plot(cor_conf[order(cor_conf)],pch=19, 
     main="correlation  with log rate cvd", ylab="cor", xlab=" ", xaxt="n")
axis(1, at=1:(length(names_conf)), labels=names_conf[order(cor_conf)], cex.axis=0.5, las=2)

# zoom in tails
par(mfrow=c(1,2))
plot(cor_conf[order(cor_conf)][1:5], pch=19, 
     main="cor.  with log rate cvd - zoom in left tail", ylab="cor", xlab=" ", xaxt="n")
axis(1, at=1:5, labels=names_conf[order(cor_conf)][1:5], cex.axis=0.6, las=2)

plot(cor_conf[order(cor_conf)][19:23],pch=19, 
     main="cor.  with log rate cvd - zoom in right tail", ylab="cor", xlab=" ", xaxt="n")
axis(1, at=1:5, labels=names_conf[order(cor_conf)][19:23], cex.axis=0.6, las=2)


################################################################################

# choosing the income
#  logMedHInc00_13

cor(data_clean$logMedHInc00_13,data_clean$Exposure)
cor(data_clean$logMedHInc00_13,data_clean$lograte)

corr_income=sapply(conf_var[-23], function(i) cor(data_clean[,i],data_clean$logMedHInc00_13))
names_proxies=names_conf[-23]

par(mfrow=c(1,1))
plot(corr_income[order(corr_income)],pch=19, 
     main="correlation with income", ylab="cor", xlab=" ", xaxt="n")
axis(1, at=1:(length(names_proxies)), labels=names_proxies[order(corr_income)], cex.axis=0.6, las=2)


# delating the varible with corr <0.2 & >-0.2
names_proxies=names_proxies[corr_income>0.2 | corr_income<(-0.2)]
corr_income=corr_income[corr_income>0.2 | corr_income<(-0.2)]

# remove variable for county level (we don't have zip code level information)
# --- "BMI2013"
names_proxies=names_proxies[-5]
corr_income=corr_income[-5]


# ----  TESTING ASSUMPIONS ----

# updated code by Kate 
con_ind_treat = sapply(c("PctWhite", "PctHighSchool", "logHVal00_13"), function(i)
  CondIndTest(data_clean[, which(colnames(data_clean) == i)],
              data_clean$Exposure,
              data_clean$logMedHInc00_13)$pvalue)

# old code
random_units=sample(1:(dim(data_clean)[1]),400)
con_ind_treat=sapply(names_proxies, function(i) 
  CondIndTest(data_clean[random_units,which(colnames(data_clean)==i)], 
              data_clean$Exposure[random_units], 
              data_clean$logMedHInc00_13[random_units])$pvalue)

con_ind_outcome=sapply(names_proxies, function(i) 
  CondIndTest(data_clean[random_units,which(colnames(data_clean)==i)], 
              data_clean$lograte[random_units], 
              data_clean$logMedHInc00_13[random_units])$pvalue)

round(con_ind_treat,4)
round(con_ind_outcome,4)

names(con_ind_treat[con_ind_treat<0.1])
names(con_ind_treat[con_ind_outcome<0.1])


##########################################################

# possible proxies

proxies=intersect(names(con_ind_treat[con_ind_treat<0.05]),names(con_ind_treat[con_ind_outcome<0.05]))

con_ind_proxy=sapply(proxies, function(i) 
  sapply(proxies[-c(1:which(proxies==i))], function(j) 
    CondIndTest(data_clean[random_units,which(colnames(data_clean)==i)], 
              data_clean[random_units,which(colnames(data_clean)==j)], 
              data_clean$logMedHInc00_13[random_units])$pvalue))


# keep PctWhite  and  PctHighSchool


dataset_final=data_clean[,c(1:3,20,24,28,4,6)]
summary(dataset_final)


plot(dataset_final[,6], dataset_final[,7])
plot(dataset_final[,6], dataset_final[,8])

plot(dataset_final[,8],dataset_final[,4])


#############################################################
