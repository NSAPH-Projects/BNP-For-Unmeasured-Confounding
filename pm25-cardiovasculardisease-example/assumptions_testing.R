# Install and load the 'dagitty' package
install.packages("dagitty")
library(dagitty)


data_clean$logPctblPvt2013 <- log(data_clean$PctblPvt2013)
data_clean$logPctblPvt2013[data_clean$logPctblPvt2013 == "-Inf"] = -4
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

                PctWhite -> logMedHInc00_13
                PctWhite -> PctblPvt2013
                PctHighSchool -> logMedHInc00_13
                PctHighSchool ->  PctblPvt2013
                logMedHInc00_13 -> logHVal00_13
                logMedHInc00_13  -> PctblPvt2013
                PctblPvt2013 -> Exposure
                PctblPvt2013 -> lograte

                PctHighSchool ->PctOwnHs2013
                logMedHInc00_13 -> PctOwnHs2013
                  PctblPvt2013 -> PctOwnHs2013
                logHVal00_13 -> PctOwnHs2013
               }"
)


dag2 <- dagitty(
  "dag {

                Exposure -> lograte
                PctWhite -> logMedHInc00_13
                PctWhite -> PctblPvt2013
                PctWhite ->logHVal00_13
                PctHighSchool -> logMedHInc00_13
                 PctHighSchool ->  PctblPvt2013
                logMedHInc00_13 -> logHVal00_13
                logMedHInc00_13  -> PctblPvt2013
                logHVal00_13 -> Exposure
                logHVal00_13 -> lograte
                PctblPvt2013 ->logHVal00_13
                logMedHInc00_13 -> logHVal00_13
                PctHighSchool ->PctOwnHs2013
                logMedHInc00_13 -> PctOwnHs2013
                  PctblPvt2013 -> PctOwnHs2013
                logHVal00_13 -> PctOwnHs2013
               }"
)



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
print(result)

ciTest("PctWhite",
       "PctblPvt2013",
       "logMedHInc00_13",
       data_clean2,
       type = "cis.loess")
