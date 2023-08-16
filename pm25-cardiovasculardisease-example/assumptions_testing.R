# Install and load the 'dagitty' package
library(dagitty)

load(
  "/n/dominici_nsaph_l3/Lab/projects/pm25-cardiovasculardisease-bnp/data_allariables_cleaned.RData"
)

data_clean$logPctblPvt2013 <- log(data_clean$PctblPvt2013 +1e10-8)
data_clean2 <- data_clean[, c(
  "Exposure",
  "PctHighSchool",
  "PctWhite",
  "logHVal00_13",
  "logMedHInc00_13",
  "lograte",
  "PctblPvt2013",
  "PctblSch2013",
  "PctOwnHs2013",
  "PctOccupied",
  'PctPoor',
  "logPctblPvt2013"
)]

# Create a DAG

dag <- dagitty(
  "dag {

                Exposure -> lograte
                logMedHInc00_13 -> Exposure
                logMedHInc00_13 -> lograte

                PctWhite -> logMedHInc00_13
                PctHighSchool -> logMedHInc00_13
                PctblSch2013 -> logMedHInc00_13

                logMedHInc00_13 -> logHVal00_13
                logMedHInc00_13  -> PctOwnHs2013

                logMedHInc00_13 -> logPctblPvt2013
                logMedHInc00_13 -> PctPoor
                logMedHInc00_13 -> PctOccupied

               }"
)


impliedConditionalIndependencies(dag)
res <- localTests(
  dag,
  data_clean,
  type = "cis.loess",
  R = 1000,
  max.conditioning.variables = 2
)
res <- res[order(abs(res$estimate)), ]
plotLocalTestResults(tail(res, 30))
print(res)

# check assumptions

# A3 about Z
  round(ciTest("PctOccupied",
               "lograte",
               c("logMedHInc00_13","Exposure"),
               data_clean2,
               type = "cis.loess", 
               R = 5000), 5)
# A4 about W
round(ciTest("PctOwnHs2013",
             "Exposure",
             c("logMedHInc00_13"),
             data_clean2,
             type = "cis.loess", 
             R = 5000), 5)
# A5
  round(ciTest("PctOwnHs2013",
       "PctOccupied",
       "logMedHInc00_13",
       data_clean2,
       type = "cis.loess", 
       R = 5000),5)

# A6 
  round(ciTest("PctOwnHs2013",
               "logMedHInc00_13",
               NULL,
               data_clean2,
               type = "cis", 
               R = 5000),5)
  
# A7 
  round(ciTest("PctOccupied",
               "logMedHInc00_13",
               "Exposure",
               data_clean2,
               type = "cis", 
               R = 5000),5)
  


