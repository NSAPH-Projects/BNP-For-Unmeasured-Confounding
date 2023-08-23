################################################################################
# ----    Check  the conditional independence assumptions A3-A5 and the 
# ----    association assumptions A6 -A7 required by our BNP-NC method 
################################################################################



# Install and load the 'dagitty' package
library(dagitty)

load(
  "/n/dominici_nsaph_l3/Lab/projects/pm25-cardiovasculardisease-bnp/data_allariables_cleaned.RData"
)

#log transform variable PctblPvt2013 because the original scale does not satisfy 
# our normality assumption on variables
data_clean$logPctblPvt2013 <- log(data_clean$PctblPvt2013 +1e10-8)
data_clean <- data_clean[, c(
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

# check assumptions A3-A7

# A3 about Z
  round(ciTest("PctOccupied",
               "lograte",
               c("logMedHInc00_13","Exposure"),
               data_clean,
               type = "cis.loess", 
               R = 5000), 5)
# A4 about W
round(ciTest("PctOwnHs2013",
             "Exposure",
             c("logMedHInc00_13"),
             data_clean,
             type = "cis.loess", 
             R = 5000), 5)
# A5
  round(ciTest("PctOwnHs2013",
       "PctOccupied",
       "logMedHInc00_13",
       data_clean,
       type = "cis.loess", 
       R = 5000),5)

# A6 
  round(ciTest("PctOwnHs2013",
               "logMedHInc00_13",
               NULL,
               data_clean,
               type = "cis", 
               R = 5000),5)
  
# A7 
  round(ciTest("PctOccupied",
               "logMedHInc00_13",
               "Exposure",
               data_clean,
               type = "cis", 
               R = 5000),5)
  