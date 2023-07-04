ReformData <- function(dataset) {
  
  dataset[, LowPov2013 := PctblPvt2013 <= 5]
  dataset[, HighPov2013 := PctblPvt2013 > 15]
  dataset[, PctblPvt2013 := NULL]
  
  # Categorizing PctOccupied
  dataset[, LowOccupied := PctOccupied <= 0.9]
  dataset[, HighOccupied := PctOccupied > 0.95]
  dataset[, PctOccupied := NULL]
  
  # Categorizing PctHisp
  dataset[, HighHispanic := PctHisp > 0.2]
  dataset[, LowHispanic := PctHisp <= 0.02]
  
  # Making covariates mean 0, and scale the continuous ones.
  dataset[, PctWhite := scale(PctWhite)]
  dataset[, PctHisp := scale(PctHisp)]
  dataset[, PctHighSchool := scale(PctHighSchool)]
  dataset[, PctPoor := scale(PctPoor)]
  dataset[, PctFemale := scale(PctFemale)]
  dataset[, PctMovedIn5 := scale(PctMovedIn5)]
  dataset[, AvgCommute := scale(AvgCommute)]
  dataset[, BMI2013 := scale(BMI2013)]
  dataset[, PctblSch2013 := scale(PctblSch2013)]
  dataset[, PctOwnHs2013 := scale(PctOwnHs2013)]
  dataset[, smokerate2013 := scale(smokerate2013)]
  dataset[, avgASOStemp := scale(avgASOStemp)]
  dataset[, avgASOSdew := scale(avgASOSdew)]
  dataset[, avgASOSrelh := scale(avgASOSrelh)]
  dataset[, mean_age := scale(mean_age)]
  dataset[, Female_rate := scale(Female_rate)]
  dataset[, Dual_rate := scale(Dual_rate)]
  dataset[, logPopPerSQM := scale(logPopPerSQM)]
  dataset[, logTotPop := scale(logTotPop)]
  dataset[, logHVal00_13 := scale(logHVal00_13)]
  dataset[, logMedHInc00_13 := scale(logMedHInc00_13)]
  dataset[, LowPov2013 := scale(LowPov2013, scale = FALSE)]
  dataset[, HighPov2013 := scale(HighPov2013, scale = FALSE)]
  dataset[, LowOccupied := scale(LowOccupied, scale = FALSE)]
  dataset[, HighOccupied := scale(HighOccupied, scale = FALSE)]
  dataset[, HighHispanic := scale(HighHispanic, scale = FALSE)]
  dataset[, LowHispanic := scale(LowHispanic, scale = FALSE)]

  return(dataset)
}