# Functions to estimate DDI effects using both a multiplicative model and an
# additive model in an SCCS design. Assumes the SCCS interval data has already
# been generated.

fitModelWithReri <- function(sccsIntervalData) {
  cd <- convertToCyclopsData(outcomes = sccsIntervalData$outcomes,
                             covariates = sccsIntervalData$covariates,
                             modelType = "cpr",
                             addIntercept = FALSE)
  fit <- fitCyclopsModel(cd, control = createControl(threads = 10))
  reri <- computeReri(fit, 1002, 1003, 1004)
  covariateIds <- c(99, 1000, 1001, 1002, 1003, 1004)
  coefs <- coef(fit)[as.character(covariateIds)]
  # If we compute all CIS in one call to confint, and one CI is invalid, we
  # will get NA for all, so calling confint once per covariate:
  cis <- list()
  for (i in seq_along(covariateIds)) {
    cis[[i]] <- confint(fit, covariateIds[i])
  }
  cis <- do.call(rbind, cis)
  estimates <- sccsIntervalData$covariateRef |>
    select(covariateId, covariateName) |>
    collect() |>
    inner_join(
      tibble(
        covariateId = covariateIds,
        irr = exp(coefs),
        ci95lb = exp(cis[, 2]),
        ci95ub = exp(cis[, 3])
      ) ,
      by = join_by("covariateId")
    )
  return(list(estimates, reri))
}

computeReri <- function(fit, covariateA = 1002, covariateB = 1003, covariateAb = 1004) {
  # Compute covariance matrix like this until issue is solved: https://github.com/OHDSI/Cyclops/issues/80
  covariates <- c(covariateA, covariateB, covariateAb)
  covariates <- Cyclops:::.checkCovariates(fit$cyclopsData, covariates)
  fi <- Cyclops:::.cyclopsGetFisherInformation(fit$cyclopsData$cyclopsInterfacePtr, covariates)
  cov_matrix <- solve(fi)

  betaA <- coef(fit)[as.character(covariateA)]
  betaB <- coef(fit)[as.character(covariateB)]
  betaAb <- coef(fit)[as.character(covariateAb)]

  # Exponentiate coefficients to compute IRRs
  irrA <- exp(betaA)
  irrB <- exp(betaB)
  irrAb <- exp(betaA + betaB + betaAb)

  # Compute Reri
  Reri <- irrAb - irrA - irrB + 1

  # Partial derivatives
  dReriDa <- exp(betaA + betaB + betaAb) - exp(betaA)
  dReriDb <- exp(betaA + betaB + betaAb) - exp(betaB)
  dReriDaB <- exp(betaA + betaB + betaAb)

  # Variance of Reri (Delta Method)
  varReri <- (
    dReriDa^2 * cov_matrix[1, 1] +
      dReriDb^2 * cov_matrix[2, 2] +
      dReriDaB^2 * cov_matrix[3, 3] +
      2 * dReriDa * dReriDb * cov_matrix[1, 2] +
      2 * dReriDa * dReriDaB * cov_matrix[1, 3] +
      2 * dReriDb * dReriDaB * cov_matrix[2, 3]
  )

  # Standard error and confidence interval
  seReri <- sqrt(varReri)
  ciReri <- c(Reri - 1.96 * seReri, Reri + 1.96 * seReri)

  result <- tibble(Reri = Reri, ci95Lb = ciReri[1], ci95Ub = ciReri[2], seReri = seReri)
  return(result)
}
