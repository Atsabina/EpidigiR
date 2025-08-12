#' Disease Prevalence Data by Region and Age Group
#'
#' A dataset containing disease prevalence data across different regions and age groups, including spatial coordinates.
#'
#' @format A data frame with 12 rows and 7 columns:
#' \describe{
#'   \item{region}{Character, region name (e.g., North, South, East, West).}
#'   \item{age_group}{Character, age group (e.g., 0-19, 20-59, 60+).}
#'   \item{cases}{Numeric, number of disease cases.}
#'   \item{population}{Numeric, population size in the region and age group.}
#'   \item{prevalence}{Numeric, prevalence percentage (cases / population * 100).}
#'   \item{lat}{Numeric, latitude for spatial mapping.}
#'   \item{lon}{Numeric, longitude for spatial mapping.}
#' }
#' @source Simulated data for demonstration purposes.
#' @examples
#' data(epi_prevalence)
#' library(sp)
#' coordinates(epi_prevalence) <- ~lon+lat
#' epi_visualize(epi_prevalence, x = "prevalence", type = "map")
#' epi_analyze(epi_prevalence, outcome = "cases", type = "summary")
"epi_prevalence"

#' SIR Model Simulation Data
#'
#' A dataset containing simulated SIR model outputs for a population of 1000.
#'
#' @format A data frame with 50 rows and 4 columns:
#' \describe{
#'   \item{time}{Numeric, time point (1 to 50 days).}
#'   \item{Susceptible}{Numeric, number of susceptible individuals.}
#'   \item{Infected}{Numeric, number of infected individuals.}
#'   \item{Recovered}{Numeric, number of recovered individuals.}
#' }
#' @source Generated using \code{epi_analyze(type = "sir", N = 1000, beta = 0.3, gamma = 0.1, days = 50)}.
#' @examples
#' data(sir_data)
#' epi_visualize(sir_data, x = "time", y = "Infected", type = "curve")
"sir_data"

#' Genomic SNP-Case Data
#'
#' A dataset containing simulated genotypes and case-control status for SNP association analysis.
#'
#' @format A data frame with 100 rows and 2 columns:
#' \describe{
#'   \item{genotypes}{Numeric, genotype (0 = AA, 1 = Aa, 2 = aa).}
#'   \item{cases}{Numeric, case (1) or control (0) status.}
#' }
#' @source Simulated data for demonstration purposes.
#' @examples
#' data(geno_data)
#' epi_model(geno_data, type = "snp")
"geno_data"

#' Machine Learning Data for Disease Risk Prediction
#'
#' A dataset containing simulated patient data for predicting disease risk, suitable for logistic regression, clustering, Random Forest, and SVM.
#'
#' @format A data frame with 100 rows and 5 columns:
#' \describe{
#'   \item{outcome}{Numeric, binary disease status (0 = healthy, 1 = diseased).}
#'   \item{age}{Numeric, patient age (years).}
#'   \item{exposure}{Numeric, exposure level (0 to 1, e.g., environmental risk).}
#'   \item{genetic_risk}{Numeric, genetic risk score (0 to 1).}
#'   \item{region}{Character, region name (e.g., North, South, East, West).}
#' }
#' @source Simulated data for demonstration purposes.
#' @examples
#' data(ml_data)
#' ml_data$outcome <- as.factor(ml_data$outcome)
#' epi_model(ml_data, formula = outcome ~ age + exposure + genetic_risk, type = "logistic")
#' epi_model(ml_data, formula = outcome ~ age + exposure + genetic_risk, type = "rf")
#' epi_visualize(ml_data, x = "age", y = "outcome", type = "scatter")
"ml_data"

#' NLP Data for Epidemiological Text Analysis
#'
#' A dataset containing simulated epidemiological text data, such as outbreak reports or health alerts, for NLP analysis.
#'
#' @format A data frame with 100 rows and 2 columns:
#' \describe{
#'   \item{id}{Character, unique identifier for each text entry.}
#'   \item{text}{Character, text content (e.g., outbreak descriptions, health reports).}
#' }
#' @source Simulated data for demonstration purposes.
#' @examples
#' data(nlp_data)
#' epi_analyze(nlp_data, outcome = NULL, type = "nlp", n = 5)
"nlp_data"

#' Clinical Trials Data for Epidemiological Analysis
#'
#' A dataset containing simulated clinical trial data for analyzing treatment outcomes, suitable for power calculations, logistic regression, Random Forest, and SVM.
#'
#' @format A data frame with 200 rows and 6 columns:
#' \describe{
#'   \item{trial_id}{Character, unique identifier for each trial participant.}
#'   \item{arm}{Character, treatment arm (e.g., Treatment, Control).}
#'   \item{outcome}{Numeric, binary outcome (0 = no response, 1 = response).}
#'   \item{age}{Numeric, patient age (years).}
#'   \item{health_score}{Numeric, baseline health score (0 to 100).}
#'   \item{dose}{Numeric, treatment dose level (e.g., 0 for control, 1 for low dose, 2 for high dose).}
#' }
#' @source Simulated data for demonstration purposes.
#' @examples
#' data(clinical_data)
#' clinical_data$outcome <- as.factor(clinical_data$outcome)
#' epi_model(clinical_data, formula = outcome ~ age + health_score + dose, type = "logistic")
#' epi_model(clinical_data, formula = outcome ~ age + health_score + dose, type = "rf")
#' epi_visualize(clinical_data, x = "arm", y = "outcome", type = "boxplot")
"clinical_data"

#' DALY Data for Global Health Burden
#'
#' A dataset containing simulated data for calculating Disability-Adjusted Life Years (DALY) in epidemiological studies.
#'
#' @format A data frame with 20 rows and 3 columns:
#' \describe{
#'   \item{group}{Character, population group (e.g., region or age group).}
#'   \item{yll}{Numeric, years of life lost due to premature mortality.}
#'   \item{yld}{Numeric, years lived with disability.}
#' }
#' @source Simulated data for demonstration purposes.
#' @examples
#' data(daly_data)
#' epi_analyze(daly_data, outcome = NULL, type = "daly")
"daly_data"

#' Survey Data for Age Standardization
#'
#' A dataset containing simulated survey data for age standardization in epidemiological studies.
#'
#' @format A data frame with 20 rows and 3 columns:
#' \describe{
#'   \item{age_group}{Character, age group (e.g., 0-19, 20-39, 40-59, 60+).}
#'   \item{rates}{Numeric, disease rates (e.g., cases per 1000).}
#'   \item{pop_weights}{Numeric, population weights for standardization.}
#' }
#' @source Simulated data for demonstration purposes.
#' @examples
#' data(survey_data)
#' epi_analyze(survey_data, outcome = NULL, type = "age_standardize")
"survey_data"

#' Diagnostic Test Data for Evaluation
#'
#' A dataset containing simulated data for evaluating diagnostic tests in epidemiological studies.
#'
#' @format A data frame with 10 rows and 5 columns:
#' \describe{
#'   \item{test_id}{Character, unique identifier for each test.}
#'   \item{true_positives}{Numeric, number of true positive results.}
#'   \item{false_positives}{Numeric, number of false positive results.}
#'   \item{true_negatives}{Numeric, number of true negative results.}
#'   \item{false_negatives}{Numeric, number of false negative results.}
#' }
#' @source Simulated data for demonstration purposes.
#' @examples
#' data(diagnostic_data)
#' epi_analyze(diagnostic_data, outcome = NULL, type = "diagnostic")
"diagnostic_data"

#' Survival Analysis Data
#'
#' A dataset containing simulated data for survival analysis in epidemiological studies.
#'
#' @format A data frame with 100 rows and 3 columns:
#' \describe{
#'   \item{id}{Character, unique identifier for each individual.}
#'   \item{time}{Numeric, time to event (e.g., years until death or censoring).}
#'   \item{status}{Numeric, event status (0 = censored, 1 = event occurred).}
#' }
#' @source Simulated data for demonstration purposes.
#' @examples
#' data(survival_data)
#' epi_model(survival_data, type = "survival")
#' epi_visualize(survival_data, x = "time", y = "status", type = "scatter")
"survival_data"
