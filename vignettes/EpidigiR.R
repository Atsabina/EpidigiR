## ----setup, include=FALSE-----------------------------------------------------
    knitr::opts_chunk$set(echo = TRUE)

# Install required dependencies if not already installed
required_packages <- c("deSolve", "sp", "tm", "glmnet", "caret", "kernlab", "survival")
for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}
# Load required packages
library(deSolve)
library(sp)
library(tm)
library(glmnet)
library(caret)
library(kernlab)
library(survival)
library(randomForest)
library(EpidigiR)
# Ensure ml_data and clinical_data outcomes are factors for caret models
data(ml_data)
ml_data$outcome <- as.factor(ml_data$outcome)
data(clinical_data)
clinical_data$outcome <- as.factor(clinical_data$outcome)

## -----------------------------------------------------------------------------
data(epi_prevalence)
epi_analyze(epi_prevalence, outcome = "cases", group = "region", type = "summary")

## -----------------------------------------------------------------------------
sir_result <- epi_analyze(data = NULL, outcome = NULL, type = "sir", N = 1000, beta = 0.3, gamma = 0.1, days = 50)
    epi_visualize(sir_result, x = "time", y = "Infected", type = "curve", main = "Epidemic Curve")

## -----------------------------------------------------------------------------
data(epi_prevalence)
coordinates(epi_prevalence) <- ~lon + lat

# Pass the Spatial object directly
epi_visualize(epi_prevalence, x = "prevalence", type = "map", main = "Prevalence Map")


## -----------------------------------------------------------------------------
data(clinical_data)
    clinical_data$outcome <- as.factor(clinical_data$outcome)
    model <- epi_model(clinical_data, formula = outcome ~ age + health_score + dose, type = "logistic")
    head(model$predictions)

## -----------------------------------------------------------------------------
rf_model <- epi_model(clinical_data, formula = outcome ~ age + health_score + dose, type = "rf")
    head(rf_model$predictions)


## -----------------------------------------------------------------------------
data(clinical_data)
epi_model(clinical_data, formula = outcome ~ age + health_score + dose, type = "logistic")

## -----------------------------------------------------------------------------
data(daly_data)
epi_analyze(daly_data, outcome = NULL, type = "daly")

