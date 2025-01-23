# Main code to run evaluation

library(SelfControlledCaseSeries)
library(readr)
source("CohortCreationFunctions.R")
source("DdiFunctions.R")
source("SetConnectionDetails.R")

exampleDdis <- read_csv("ExampleDdis.csv", show_col_types = FALSE)

drugAcohortId <- 1
drugBcohortId <- 2
drugAandBcohortId <- 100
drugAandOtherDrugsInClasOfbCohortId <-101
drugBandOtherDrugsInClasOfaCohortId <- 102

resultRows <- list()
dbi = 5
i = 1
for (dbi in 1:nrow(databases)) {
  database <- databases[dbi, ]
  writeLines(sprintf("*** Evaluating in %s ***", database$name))
  if (!file.exists(database$folder))
    dir.create(database$folder, recursive = TRUE)
  for (i in 1:nrow(exampleDdis)) {
    exampleDdi <- exampleDdis[i, ]
    writeLines(sprintf("** Evaluating example %d: %s and %s for %s **",
                       i,
                       exampleDdi$drugAname,
                       exampleDdi$drugBname,
                       exampleDdi$outcomeName))
    idString <- sprintf("_a%d_b%d_o%d", exampleDdi$drugAconceptId, exampleDdi$drugBconceptId, exampleDdi$outcomeCohortId)
    resultRowFileName <- file.path(database$folder, sprintf("resultRow%s.rds", idString))
    if (file.exists(resultRowFileName)) {
      resultRow <- readRDS(resultRowFileName)
    } else {
      otherDrugsInClassOfaConceptIds <- as.numeric(strsplit(exampleDdi$otherDrugsInClassOfaConceptIds, ";")[[1]])
      otherDrugsInClassOfbConceptIds <- as.numeric(strsplit(exampleDdi$otherDrugsInClassOfbConceptIds, ";")[[1]])

      createCohorts(
        connectionDetails,
        database$cdmDatabaseSchema,
        database$cohortDatabaseSchema,
        database$cohortTable,
        drugAconceptId = exampleDdi$drugAconceptId,
        drugBconceptId = exampleDdi$drugBconceptId,
        otherDrugsInClassOfaConceptIds = otherDrugsInClassOfaConceptIds,
        otherDrugsInClassOfbConceptIds = otherDrugsInClassOfbConceptIds,
        outcomeCohortId = exampleDdi$outcomeCohortId,
        drugAcohortId = drugAcohortId,
        drugBcohortId = drugBcohortId,
        drugAandBcohortId = drugAandBcohortId,
        drugAandOtherDrugsInClasOfbCohortId = drugAandOtherDrugsInClasOfbCohortId,
        drugBandOtherDrugsInClasOfaCohortId = drugBandOtherDrugsInClasOfaCohortId,
        surveillanceDays = 14
      )

      sccsData <- getDbSccsData(connectionDetails = connectionDetails,
                                cdmDatabaseSchema = database$cdmDatabaseSchema,
                                outcomeDatabaseSchema = database$cohortDatabaseSchema,
                                outcomeTable = database$cohortTable,
                                outcomeIds = exampleDdi$outcomeCohortId,
                                exposureDatabaseSchema = database$cohortDatabaseSchema,
                                exposureTable = database$cohortTable,
                                exposureIds = c(drugAcohortId,
                                                drugBcohortId,
                                                drugAandBcohortId,
                                                drugAandOtherDrugsInClasOfbCohortId,
                                                drugBandOtherDrugsInClasOfaCohortId),
                                studyStartDates = "20140101",
                                studyEndDates = "20241231",
                                maxCasesPerOutcome = 100000)
      saveSccsData(sccsData, file.path(database$folder, sprintf("sccsData%s.zip", idString)))
      sccsData <- loadSccsData(file.path(database$folder, sprintf("sccsData%s.zip", idString)))

      studyPop <- createStudyPopulation(sccsData = sccsData,
                                        outcomeId = exampleDdi$outcomeCohortId,
                                        firstOutcomeOnly = TRUE,
                                        naivePeriod = 180)
      preExposure1 <- createEraCovariateSettings(label = "Pre drug A exposure",
                                                 includeEraIds = drugAcohortId,
                                                 start = -30,
                                                 end = -1,
                                                 endAnchor = "era start")
      preExposure2 <- createEraCovariateSettings(label = "Pre drug B exposure",
                                                 includeEraIds = drugBcohortId,
                                                 start = -30,
                                                 end = -1,
                                                 endAnchor = "era start")
      covar1 <- createEraCovariateSettings(label = "Drug A",
                                           includeEraIds = drugAcohortId,
                                           start = 1,
                                           end = 0,
                                           endAnchor = "era end")
      covar2 <- createEraCovariateSettings(label = "Drug B",
                                           includeEraIds = drugBcohortId,
                                           start = 1,
                                           end = 0,
                                           endAnchor = "era end")
      covarBoth <- createEraCovariateSettings(label = "Both",
                                              includeEraIds = drugAandBcohortId,
                                              start = 1,
                                              end = 0,
                                              endAnchor = "era end")
      covar1AndOther<- createEraCovariateSettings(label = "Drug A and other drugs in the class of B",
                                                  includeEraIds = drugAandOtherDrugsInClasOfbCohortId,
                                                  start = 0,
                                                  end = 0,
                                                  endAnchor = "era end")
      covar2AndOther<- createEraCovariateSettings(label = "Drug B and other drugs in the class of A",
                                                  includeEraIds = drugBandOtherDrugsInClasOfaCohortId,
                                                  start = 0,
                                                  end = 0,
                                                  endAnchor = "era end")
      seasonalityCovariateSettings <- createSeasonalityCovariateSettings()
      calendarTimeSettings <- createCalendarTimeCovariateSettings()
      sccsIntervalData <- createSccsIntervalData(studyPopulation = studyPop,
                                                 sccsData = sccsData,
                                                 eraCovariateSettings = list(preExposure1,
                                                                             preExposure2,
                                                                             covar1,
                                                                             covar2,
                                                                             covarBoth,
                                                                             covar1AndOther,
                                                                             covar2AndOther),
                                                 seasonalityCovariateSettings = seasonalityCovariateSettings,
                                                 calendarTimeCovariateSettings = calendarTimeSettings)
      # attr(sccsIntervalData, "metaData")$covariateStatistics

      estimates <- fitModelWithReri(sccsIntervalData)

      estimateA <- estimates[[1]] |>
        filter(covariateId == 1002) |>
        transmute(irrA = sprintf("%0.2f (%0.2f-%0.2f)", irr, ci95lb, ci95ub))
      estimateB <- estimates[[1]] |>
        filter(covariateId == 1003) |>
        transmute(irrB = sprintf("%0.2f (%0.2f-%0.2f)", irr, ci95lb, ci95ub))
      estimateAandB <- estimates[[1]] |>
        filter(covariateId == 1004) |>
        transmute(irrAandB = sprintf("%0.2f (%0.2f-%0.2f)", irr, ci95lb, ci95ub))
      reri <- estimates[[2]] |>
        transmute(reriAandB = sprintf("%0.2f (%0.2f-%0.2f)", Reri, ci95Lb, ci95Ub))

      # Compute SCCS diagnostics
      if (!1000 %in% estimates[[1]]$covariateId) {
        preExposure <- tibble(irr = NA, ci95lb = NA, ci95ub = NA)
      } else {
        preExposure <- estimates[[1]] |>
          filter(covariateId == 1000) |>
          select(irr,  ci95lb, ci95ub)
      }
      preExposureA <- preExposure |>
        transmute(irrPreExposureA = sprintf("%0.2f (%0.2f-%0.2f)", irr, ci95lb, ci95ub),
                  passPreExposureA = (is.na(ci95lb) | is.na(ci95ub)) || (ci95lb < 1.25 && ci95ub > 0.8))

      if (!1001 %in% estimates[[1]]$covariateId) {
        preExposure <- tibble(irr = NA, ci95lb = NA, ci95ub = NA)
      } else {
        preExposure <- estimates[[1]] |>
          filter(covariateId == 1001) |>
          select(irr,  ci95lb, ci95ub)
      }
      preExposureB <- preExposure |>
        transmute(irrPreExposureB = sprintf("%0.2f (%0.2f-%0.2f)", irr, ci95lb, ci95ub),
                  passPreExposureB = (is.na(ci95lb) | is.na(ci95ub)) || (ci95lb < 1.25 && ci95ub > 0.8))

      if (!99 %in% estimates[[1]]$covariateId) {
        endOfObs <- tibble(irr = NA, ci95lb = NA, ci95ub = NA)
      } else {
        endOfObs <- estimates[[1]] |>
          filter(covariateId == 99) |>
          select(irr,  ci95lb, ci95ub)
      }
      endOfObs <- endOfObs |>
        transmute(irrEndOfObservation = sprintf("%0.2f (%0.2f-%0.2f)", irr, ci95lb, ci95ub),
                  passEndOfObservation = (is.na(ci95lb) | is.na(ci95ub)) || (ci95lb < 2 && ci95ub > 0.5))

      # Need proper SccsModel object to compute time trend, so fitting again:
      model <- fitSccsModel(sccsIntervalData, control = createControl(threads = 10), profileBounds = NULL)
      timeStability <- computeTimeStability(studyPopulation = studyPop, sccsModel = model) |>
        select(pTimeStability = p, passTimeStability = stable)

      rareOutcome <- checkRareOutcomeAssumption(studyPopulation = studyPop) |>
        mutate(rare = if_else(is.na(rare), TRUE, rare)) |>
        select(outcomeProportion, passRareOutcome = rare)

      resultRow <- exampleDdi |>
        select(drugAEffect, drugBeffect, interaction, drugAname, drugBname, outcomeName) |>
        mutate(database = database$name) |>
        bind_cols(estimateA,
                  estimateB,
                  estimateAandB,
                  reri,
                  preExposureA,
                  preExposureB,
                  endOfObs,
                  timeStability,
                  rareOutcome)
      saveRDS(resultRow, resultRowFileName)
    }
    resultRows[[length(resultRows) + 1]] <- resultRow
  }
}

resultRows <- bind_rows(resultRows)
write_csv(resultRows, "EvaluationResults.csv")
