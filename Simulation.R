# Simulation to test Reri confidence interval coverage

library(Cyclops)
library(survival)

source("DdiFunctions.R")

simulateOne <- function(seed) {
  set.seed(seed)

  nPersons <- 1000
  trueIrrA <- 0.5
  trueIrrB <- 2.2
  trueIrrAB <- 3.7

  trueReri <- trueIrrAB- trueIrrA - trueIrrB + 1
  intervals <- tibble(
    personId = rep(seq(nPersons), each = 4),
    baselineRate = rep(runif(nPersons, 0.01, 0.05), each = 4),
    a = rep(c(0, 1, 0, 1), nPersons),
    b = rep(c(0, 0, 1, 1), nPersons),
    time = runif(nPersons * 4, 1, 10)
  )
  intervals <- intervals |>
    mutate(rate = baselineRate * case_when(a == 1 & b == 0 ~ trueIrrA,
                                           a == 0 & b == 1 ~ trueIrrB,
                                           a == 1 & b == 1 ~ trueIrrAB,
                                           TRUE ~ 1)) |>
    mutate(outcomes = rpois(n(), rate * time))
  sum(intervals$outcomes)

  cyclopsData <- createCyclopsData(outcomes ~ a + b + a*b + strata(personId) + offset(log(time)),
                                   data = intervals,
                                   modelType = "cpr")
  fit <- fitCyclopsModel(cyclopsData)
  Reri <- computeReri(fit, "a", "b", "a:b")
  return(trueReri > Reri$ci95Lb & trueReri < Reri$ci95Ub)
}
coverage <- sapply(1:1000, simulateOne)
mean(coverage)
