library(Capr)
library(dplyr)

createCohorts <- function(
    connectionDetails,
    cdmDatabaseSchema,
    cohortDatabaseSchema,
    cohortTable,
    drugAconceptId,
    drugBconceptId,
    otherDrugsInClassOfaConceptIds,
    otherDrugsInClassOfbConceptIds,
    outcomeCohortId,
    drugAcohortId = 1,
    drugBcohortId = 2,
    drugAandBcohortId = 100,
    drugAandOtherDrugsInClasOfbCohortId = 101,
    drugBandOtherDrugsInClasOfaCohortId = 102,
    surveillanceDays = 14) {

  otherDrugsInClassOfaCohortId <- 999
  otherDrugsInClassOfbCohortId <- 998

  # drugAconceptId <- 1310149 # Warfarin
  # drugBconceptId <- 1309944 # Amiodarone
  # otherDrugsInClassOfaConceptIds <- c(19024063, # Acenocoumarol
  #                                     19035344, # Phenprocoumon
  #                                     1325124)  # Dicoumarol
  #
  # otherDrugsInClassOfbConceptIds <- c(40163615) # Dronedarone
  # outcomecohortId <- 77 # GI bleed
  # drugAcohortId <- 1
  # drugBcohortId <- 2
  # otherDrugsInClassOfaConceptIds <- 3
  # otherDrugsInClassOfbConceptIds <- 4
  # drugAandBcohortId <- 100
  # drugAandOtherDrugsInClasOfbCohortId  <- 101
  # druBAandOtherDrugsInClasOfbCohortId  <- 102

  connection <- DatabaseConnector::connect(connectionDetails)
  on.exit(DatabaseConnector::disconnect(connection))

  message("Creating cohort tables")
  cohortTableNames <- CohortGenerator::getCohortTableNames(cohortTable)
  CohortGenerator::createCohortTables(connection = connection,
                                      cohortDatabaseSchema = cohortDatabaseSchema,
                                      cohortTableNames = cohortTableNames)

  message("Generating base cohorts")
  drugAcs <- cs(
    descendants(drugAconceptId),
    name = "Drug A"
  )
  drugAcohort <- cohort(
    entry = entry(
      drugExposure(drugAcs)
    ),
    exit = exit(
      endStrategy = drugExit(drugAcs,
                             persistenceWindow = 30,
                             surveillanceWindow = 0)
    )
  )
  drugBcs <- cs(
    descendants(drugBconceptId),
    name = "Drug B"
  )
  drugBcohort <- cohort(
    entry = entry(
      drugExposure(drugBcs)
    ),
    exit = exit(
      endStrategy = drugExit(drugBcs,
                             persistenceWindow = 30,
                             surveillanceWindow = 0)
    )
  )
  otherDrugsInClassAcs <- cs(
    descendants(otherDrugsInClassOfaConceptIds),
    name = "Other drugs in class of A"
  )
  otherDrugsInClassACohort <- cohort(
    entry = entry(
      drugExposure(otherDrugsInClassAcs)
    ),
    exit = exit(
      endStrategy = drugExit(otherDrugsInClassAcs,
                             persistenceWindow = 30,
                             surveillanceWindow = 0)
    )
  )
  otherDrugsInClassBcs <- cs(
    descendants(otherDrugsInClassOfbConceptIds),
    name = "Other drugs in class of B"
  )
  otherDrugsInClassBCohort <- cohort(
    entry = entry(
      drugExposure(otherDrugsInClassBcs)
    ),
    exit = exit(
      endStrategy = drugExit(otherDrugsInClassBcs,
                             persistenceWindow = 30,
                             surveillanceWindow = 0)
    )
  )
  cohortDefinitionSet <- makeCohortSet(drugAcohort,
                                       drugBcohort,
                                       otherDrugsInClassACohort,
                                       otherDrugsInClassBCohort)
  cohortDefinitionSet$cohortId <- c(drugAcohortId,
                                    drugBcohortId,
                                    otherDrugsInClassOfaCohortId,
                                    otherDrugsInClassOfbCohortId)
  cohortDefinitionSet <- bind_rows(cohortDefinitionSet,
                                   PhenotypeLibrary::getPlCohortDefinitionSet(outcomeCohortId))
  CohortGenerator::generateCohortSet(connection = connection,
                                     cdmDatabaseSchema = cdmDatabaseSchema,
                                     cohortDatabaseSchema = cohortDatabaseSchema,
                                     cohortTableNames = cohortTableNames,
                                     cohortDefinitionSet = cohortDefinitionSet)

  message("Generating interaction cohorts")
  sql <- "
  DELETE FROM @cohort_database_schema.@cohort_table
  WHERE cohort_definition_id IN (@cohort_1_and_2_id);

  INSERT INTO @cohort_database_schema.@cohort_table (cohort_definition_id, subject_id, cohort_start_date, cohort_end_date)
  SELECT @cohort_1_and_2_id AS cohort_definition_id,
  	cohort_1.subject_id,
  	CASE WHEN
  	  cohort_1.cohort_start_date >= cohort_2.cohort_start_date THEN cohort_1.cohort_start_date
  	  ELSE cohort_2.cohort_start_date
  	END AS cohort_start_date,
  	DATEADD(DAY, @surveillance_days,
    	CASE WHEN
    	  cohort_1.cohort_end_date <= cohort_2.cohort_end_date THEN cohort_1.cohort_end_date
    	  ELSE cohort_2.cohort_end_date
    	END
    ) AS cohort_end_date
  FROM @cohort_database_schema.@cohort_table cohort_1
  INNER JOIN  @cohort_database_schema.@cohort_table cohort_2
    ON cohort_1.subject_id = cohort_2.subject_id
      AND cohort_1.cohort_start_date <= cohort_2.cohort_end_date
      AND cohort_1.cohort_end_date >= cohort_2.cohort_start_date
  WHERE cohort_1.cohort_definition_id = @cohort_1_id
    AND cohort_2.cohort_definition_id = @cohort_2_id;
  "
  DatabaseConnector::renderTranslateExecuteSql(connection = connection,
                                               sql = sql,
                                               cohort_database_schema = cohortDatabaseSchema,
                                               cohort_table = cohortTable,
                                               cohort_1_id = drugAcohortId,
                                               cohort_2_id = drugBcohortId,
                                               cohort_1_and_2_id = drugAandBcohortId,
                                               surveillance_days = surveillanceDays)

  DatabaseConnector::renderTranslateExecuteSql(connection = connection,
                                               sql = sql,
                                               cohort_database_schema = cohortDatabaseSchema,
                                               cohort_table = cohortTable,
                                               cohort_1_id = drugAcohortId,
                                               cohort_2_id = otherDrugsInClassOfaCohortId,
                                               cohort_1_and_2_id = drugAandOtherDrugsInClasOfbCohortId,
                                               surveillance_days = 0)

  DatabaseConnector::renderTranslateExecuteSql(connection = connection,
                                               sql = sql,
                                               cohort_database_schema = cohortDatabaseSchema,
                                               cohort_table = cohortTable,
                                               cohort_1_id = drugBcohortId,
                                               cohort_2_id = otherDrugsInClassOfbCohortId,
                                               cohort_1_and_2_id = drugBandOtherDrugsInClasOfaCohortId,
                                               surveillance_days = 0)

  # Check number of subjects per cohort:
  sql <- "SELECT cohort_definition_id,
    COUNT(*) AS cohort_count,
    SUM(CAST(DATEDIFF(DAY, cohort_start_date, cohort_end_date) + 1 AS BIGINT)) AS day_count,
    MIN(cohort_start_date) AS min_date,
    MAX(cohort_start_date) AS max_date
  FROM @cohortDatabaseSchema.@cohortTable
  GROUP BY cohort_definition_id;"
  counts <- DatabaseConnector::renderTranslateQuerySql(connection = connection,
                                             sql = sql,
                                             cohortDatabaseSchema = cohortDatabaseSchema,
                                             cohortTable = cohortTable)
  invisible(counts)
}
