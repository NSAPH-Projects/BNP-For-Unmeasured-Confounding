DataClean <- function(dataset, trt = 'avgPM', max_exposure = NULL, plotCor = FALSE,
                      drop_missing = FALSE) {
  # It takes the dataset, keeps only numeric and integer covariates, calculates the
  # pairwise correlations and drops covariates of high correlation. Then it drops
  # some variables that I do not need.
  #
  # Args:
  #  dataset:   The zipcode level data
  #  trt:       The name of the PM information that is the treatment.
  #  max_exposure: If not NULL, any observations with exposure value larger than
  #  this will be dropped.
  #  plotCor:   Whether we want to plot the correlations. (FALSE)
  #  drop_missing: Drop observations with missing data.
  #
  # Returns the same dataset with fewer covariates after dropping the ones we do not
  # need.
  
  setnames(dataset, trt, 'Exposure')
  
  wh <- which(with(dataset, is.na(Exposure)))
  print(paste('Dropping', length(wh), 'ZIP codes because of missing exposure.'))
  dataset <- subset(dataset, !is.na(Exposure))
  
  if (!is.null(max_exposure)) {
    wh <- which(with(dataset, Exposure > max_exposure))
    print(paste('Dropping', length(wh), 'ZIP codes because exposure is larger than',
                max_exposure))
    dataset <- subset(dataset, Exposure <= max_exposure)
  }
  
  dataset[, Condition := NULL]
  dataset[, year := NULL]
  
  # Dropping character variables.
  dataset[, City.zip := NULL]
  dataset[, State.zip := NULL]
  dataset[, InUSE := NULL]
  dataset[, OBJECTID := NULL]
  
  # Defining the outcome.
  dataset[, lograte := log(Tot_num / Personal_years)]
  
  # Dropping the pollution variables from Yun.
  dataset[, c('PM25_two_year_mean', 'PM25_two_year_min', 'PM25_two_year_max',
              'PM25_two_year_total_point') := NULL]
  
  # Keeping only the chosen treatment.  
  PMcovs <- c('avgPM', 'meanMaxPM', 'meanMedPM', 'mean4maxPM')
  drop_PM <- setdiff(PMcovs, trt)
  dataset[, (drop_PM) := NULL]
  
  # Dropping number of monitors.
  dataset[, c('Total_monitor', 'nPMmonitors', 'nASOSstations',
              'nAQSTemp_monitors') := NULL]
  
  # Dropping ASOS altitude because areas in NH/Vt have significantly lower values.
  dataset[, avgASOSalti := NULL]
  
  # Keeping ASOS temperature info since AQS has more missing.
  dataset[, c('avgAQStemp', 'meanMaxAQStemp', 'meanMedAQStemp') := NULL]
  
  # Keeping pctUrban instead of:
  dataset[, c('PctRural', 'PctUrbanHUs', 'PctRuralHUs', 'PctInUAs') := NULL]
  dataset[, c('Rural', 'RuralHUs', 'Urban', 'UrbanHUs') := NULL]
  
  # Keeping PctHisp instead of:
  dataset[, c('HispPop', 'HspPct2013') := NULL]  # corr 0.83, 0.92
  
  # Keeping TotPop instead of:
  dataset[, c('InUAs', 'HighSchool', 'Over25', 'PovUniverse', 'Female', 'TotHUs',
              'Occupied', 'MovedIn5') := NULL]
  
  # Dropping numerator of log rate and correlated variable:
  dataset[, c('Tot_num', 'Total_den') := NULL]
  
  # Keeping PctWhite instead of
  dataset[, c('White_rate', 'BlkPct2013', 'Black_rate', 'PctBlack') := NULL]
  
  # Taking the logarithm:
  dataset[, logPopDen := log(PopDen2013)]
  dataset[, logPopPerSQM := log(PopPerSQM + 1)]
  dataset[, logLandSQM := log(LandSQMI + 1)]
  dataset[, logTotPop := log(TotPop + 1)]
  dataset[, c('PopDen2013', 'PopPerSQM', 'LandSQMI', 'TotPop') := NULL]
  
  # Reforming Household value and income
  dataset[, logHVal00_13 := log(MedianHValue + MdVlHs2013)]
  dataset[, c('MedianHValue', 'MdVlHs2013') := NULL]
  dataset[, logMedHInc00_13 := log(MedianHHInc + MdHsIcm2013)]
  dataset[, c('MedianHHInc', 'MdHsIcm2013') := NULL]
  
  # Dropping for pathway reasons:
  dataset[, LungCancerRate2013 := NULL]
  
  # Dropping the following but I dont know why.
  dataset[, c('White1', 'Black1', 'Poor', 'InUCs', 'Personal_years') := NULL]
  
  # Dropping variables because they're almost constant:
  dataset[, c('PctEye2013', 'Pctmam2013', 'PctA1c2013', 'PctAml2013',
              'PctLDL2013', 'PctUrban', 'PctInUCs') := NULL]
  
  # Keeping logPopPerSQM instead of:
  dataset[, c('logPopDen', 'logLandSQM') := NULL]
  
  # Dropping observations with missing data.
  if (drop_missing) {
    wh <- which(apply(dataset, 1, function(x) sum(is.na(x))) > 0)
    print(paste('Dropping', length(wh), 'zip codes due to missing data. This',
                'corresponds to', round(length(wh) / nrow(dataset), 4), 'of observations.'))
    dataset <- dataset[- wh]
  }
  
  if (plotCor) {
    C <- cor(dataset, use = 'pairwise.complete.obs')
    corrplot(C, tl.cex = 0.6)
  }
  
  return(dataset)
}