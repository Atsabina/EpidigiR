#' Comprehensive Epidemiological Analysis
#'
#' Performs summary statistics, SIR modeling, DALY calculation, age standardization, diagnostic test evaluation, or NLP keyword extraction.
#'
#' @param data Input data frame with relevant columns (e.g., cases, population, yll, yld, text).
#' @param outcome Outcome column name (character, e.g., "cases").
#' @param group Grouping column name (character, e.g., "region", optional).
#' @param type Analysis type: "summary", "sir", "daly", "age_standardize", "diagnostic", "nlp".
#' @param ... Additional parameters (e.g., N, beta, gamma for SIR).
#' @return A data frame with analysis results.
#' @export
epi_analyze <- function(data, outcome, group = NULL, type = c("summary", "sir", "daly", "age_standardize", "diagnostic", "nlp"), ...) {
  type <- match.arg(type)
  switch(type,
         summary = {
           if (is.null(data) || is.null(outcome)) stop("Data and outcome must be provided for summary analysis")
           if (!outcome %in% names(data)) stop("Outcome not found in data")
           if (!is.null(group) && !group %in% names(data)) stop("Group not found in data")
           stats <- if (is.null(group)) {
             data.frame(
               mean_outcome = mean(data[[outcome]], na.rm = TRUE),
               incidence_rate = sum(data[[outcome]], na.rm = TRUE) / sum(data$population, na.rm = TRUE) * 1000,
               prevalence = mean(data$prevalence, na.rm = TRUE)
             )
           } else {
             aggregate(data[[outcome]], list(data[[group]]), function(x) c(
               mean_outcome = mean(x, na.rm = TRUE),
               incidence_rate = sum(x, na.rm = TRUE) / sum(data$population[data[[group]] == unique(data[[group]])[1]], na.rm = TRUE) * 1000,
               prevalence = mean(data$prevalence[data[[group]] == unique(data[[group]])[1]], na.rm = TRUE)
             ))
           }
           stats
         },
         sir = {
           if (!requireNamespace("deSolve", quietly = TRUE)) {
             stop("Package 'deSolve' is required for SIR modeling. Please install it using install.packages('deSolve').")
           }
           params <- list(...)
           N <- if (!is.null(params$N)) params$N else 1000
           beta <- if (!is.null(params$beta)) params$beta else 0.3
           gamma <- if (!is.null(params$gamma)) params$gamma else 0.1
           days <- if (!is.null(params$days)) params$days else 50
           init <- c(S = N - 1, I = 1, R = 0)
           parameters <- c(beta = beta, gamma = gamma, N = N)
           times <- seq(1, days, by = 1)
           sir_equations <- function(time, state, parameters) {
             with(as.list(c(state, parameters)), {
               dS <- -beta * S * I / N
               dI <- beta * S * I / N - gamma * I
               dR <- gamma * I
               list(c(dS, dI, dR))
             })
           }
           out <- deSolve::ode(y = init, times = times, func = sir_equations, parms = parameters)
           data.frame(time = out[, "time"], Susceptible = out[, "S"], Infected = out[, "I"], Recovered = out[, "R"])
         },
         daly = {
           if (is.null(data)) stop("Data must be provided for DALY analysis")
           if (!all(c("yll", "yld") %in% names(data))) stop("Data must contain 'yll' and 'yld' columns")
           data.frame(group = data$group, daly = data$yll + data$yld)
         },
         age_standardize = {
           if (is.null(data)) stop("Data must be provided for age standardization")
           if (!all(c("rates", "pop_weights") %in% names(data))) stop("Data must contain 'rates' and 'pop_weights' columns")
           standardized_rate <- sum(data$rates * data$pop_weights) / sum(data$pop_weights)
           data.frame(standardized_rate = standardized_rate)
         },
         diagnostic = {
           if (is.null(data)) stop("Data must be provided for diagnostic test evaluation")
           if (!all(c("true_positives", "false_positives", "true_negatives", "false_negatives") %in% names(data))) {
             stop("Data must contain 'true_positives', 'false_positives', 'true_negatives', 'false_negatives' columns")
           }
           sensitivity <- data$true_positives / (data$true_positives + data$false_negatives)
           specificity <- data$true_negatives / (data$true_negatives + data$false_positives)
           accuracy <- (data$true_positives + data$true_negatives) / (rowSums(data[, c("true_positives", "true_negatives", "false_positives", "false_negatives")]))
           data.frame(test_id = data$test_id, sensitivity = sensitivity, specificity = specificity, accuracy = accuracy)
         },
         nlp = {
           if (!requireNamespace("tm", quietly = TRUE)) {
             stop("Package 'tm' is required for NLP analysis. Please install it using install.packages('tm').")
           }
           params <- list(...)
           n <- if (!is.null(params$n)) params$n else 10
           corpus <- tm::Corpus(tm::VectorSource(data$text))
           corpus <- tm::tm_map(corpus, tm::content_transformer(tolower))
           corpus <- tm::tm_map(corpus, tm::removePunctuation)
           corpus <- tm::tm_map(corpus, tm::removeNumbers)
           corpus <- tm::tm_map(corpus, tm::removeWords, tm::stopwords("english"))
           tdm <- tm::TermDocumentMatrix(corpus)
           m <- as.matrix(tdm)
           word_freq <- sort(rowSums(m), decreasing = TRUE)
           data.frame(word = names(word_freq)[1:n], frequency = word_freq[1:n])
         }
  )
}

#' Unified Epidemiological Modeling
#'
#' Performs clinical trial power calculation, survival analysis, SNP association, logistic regression, k-means clustering, Random Forest, or SVM.
#'
#' @param data Input data frame with relevant columns (e.g., outcome, genotypes).
#' @param formula Model formula (optional, for survival/logistic/rf/svmRadial, e.g., "outcome ~ x").
#' @param type Model type: "power", "survival", "snp", "logistic", "kmeans", "rf", "svmRadial".
#' @param ... Additional parameters (e.g., n, effect_size for power; k for kmeans).
#' @return A data frame or list with model results.
#' @export
epi_model <- function(data, formula = NULL, type = c("power", "survival", "snp", "logistic", "kmeans", "rf", "svmRadial"), ...) {
  type <- match.arg(type)
  switch(type,
         power = {
           params <- list(...)
           n <- if (!is.null(params$n)) params$n else 100
           effect_size <- if (!is.null(params$effect_size)) params$effect_size else 0.5
           sd <- if (!is.null(params$sd)) params$sd else 1
           power <- power.t.test(n = n, delta = effect_size, sd = sd, sig.level = 0.05, power = NULL)$power
           data.frame(n = n, effect_size = effect_size, sd = sd, power = power)
         },
         survival = {
           if (!requireNamespace("survival", quietly = TRUE)) {
             stop("Package 'survival' is required for survival analysis. Please install it using install.packages('survival').")
           }
           if (is.null(data) || !all(c("time", "status") %in% names(data))) stop("Data must contain 'time' and 'status' columns")
           surv_obj <- survival::Surv(data$time, data$status)
           fit <- survival::survfit(surv_obj ~ 1)
           list(survfit = fit, summary = summary(fit))
         },
         snp = {
           if (is.null(data) || !all(c("genotypes", "cases") %in% names(data))) stop("Data must contain 'genotypes' and 'cases' columns")
           tbl <- table(data$genotypes, data$cases)
           chi_test <- chisq.test(tbl)
           data.frame(statistic = chi_test$statistic, p_value = chi_test$p.value)
         },
         logistic = {
           if (!requireNamespace("glmnet", quietly = TRUE)) {
             stop("Package 'glmnet' is required for logistic regression. Please install it using install.packages('glmnet').")
           }
           if (is.null(data) || is.null(formula)) stop("Data and formula must be provided for logistic regression")
           mf <- model.frame(formula, data)
           y <- model.response(mf)
           x <- model.matrix(formula, data)[, -1]
           fit <- glmnet::cv.glmnet(x, y, family = "binomial")
           list(coefficients = coef(fit, s = "lambda.min"), predictions = predict(fit, newx = x, s = "lambda.min", type = "response"))
         },
         kmeans = {
           params <- list(...)
           k <- if (!is.null(params$k)) params$k else 3
           if (is.null(data)) stop("Data must be provided for k-means clustering")
           fit <- kmeans(data, centers = k)
           list(clusters = fit$cluster, centers = fit$centers)
         },
         rf = {
           if (!requireNamespace("caret", quietly = TRUE)) {
             stop("Package 'caret' is required for Random Forest. Please install it using install.packages('caret').")
           }
           if (is.null(data) || is.null(formula)) stop("Data and formula must be provided for Random Forest")
           # Ensure outcome is a factor for classification
           outcome_name <- all.vars(formula)[1]
           data[[outcome_name]] <- as.factor(data[[outcome_name]])
           fit <- caret::train(formula, data = data, method = "rf", trControl = caret::trainControl(method = "cv", number = 5), tuneLength = 3)
           list(model = fit, performance = try(caret::confusionMatrix(predict(fit, data), data[[outcome_name]]), silent = TRUE))
         },
         svmRadial = {
           if (!requireNamespace("caret", quietly = TRUE)) {
             stop("Package 'caret' is required for SVM. Please install it using install.packages('caret').")
           }
           if (is.null(data) || is.null(formula)) stop("Data and formula must be provided for SVM")
           # Ensure outcome is a factor for classification
           outcome_name <- all.vars(formula)[1]
           data[[outcome_name]] <- as.factor(data[[outcome_name]])
           fit <- caret::train(formula, data = data, method = "svmRadial", trControl = caret::trainControl(method = "cv", number = 5), tuneLength = 3)
           list(model = fit, performance = try(caret::confusionMatrix(predict(fit, data), data[[outcome_name]]), silent = TRUE))
         }
  )
}

#' Flexible Epidemiological Visualization
#'
#' Creates visualizations for prevalence mapping, epidemic curves, or general plots (scatter, boxplot).
#'
#' @param data Input data frame or SpatialPolygonsDataFrame with relevant columns.
#' @param x X-axis column name (character, e.g., "region").
#' @param y Y-axis column name (character, e.g., "prevalence", optional).
#' @param type Plot type: "map", "curve", "scatter", "boxplot".
#' @param ... Additional plotting parameters (e.g., main, xlab).
#' @return A plot (spplot for maps, base R for others).
#' @export
epi_visualize <- function(data, x, y = NULL, type = c("map", "curve", "scatter", "boxplot"), ...) {
  if (!requireNamespace("sp", quietly = TRUE)) {
    stop("Package 'sp' is required for visualization. Please install it using install.packages('sp').")
  }
  type <- match.arg(type)
  params <- list(...)
  # Define defaults
  default_main <- switch(type, map = "Disease Prevalence Map", curve = "Epidemic Curve", scatter = "Scatter Plot", boxplot = "Boxplot")
  default_xlab <- x
  default_ylab <- if (!is.null(y)) y else "Value"
  # Prioritize user-specified arguments, fallback to defaults
  main <- if (!is.null(params$main)) params$main else default_main
  xlab <- if (!is.null(params$xlab)) params$xlab else default_xlab
  ylab <- if (!is.null(params$ylab)) params$ylab else default_ylab
  # Filter out main, xlab, ylab from params to avoid duplicates
  params <- params[!names(params) %in% c("main", "xlab", "ylab")]
  switch(type,
         map = {
           if (!inherits(data, "Spatial")) stop("Data must be a Spatial object for map visualization")
           do.call(sp::spplot, c(list(obj = data, zcol = x, main = main, xlab = xlab, ylab = ylab), params))
         },
         curve = {
           if (is.null(y)) stop("Y variable must be provided for curve visualization")
           do.call(plot, c(list(x = data[[x]], y = data[[y]], type = "l", main = main, xlab = xlab, ylab = ylab), params))
         },
         scatter = {
           if (is.null(y)) stop("Y variable must be provided for scatter visualization")
           do.call(plot, c(list(x = data[[x]], y = data[[y]], main = main, xlab = xlab, ylab = ylab), params))
         },
         boxplot = {
           if (is.null(y)) stop("Y variable must be provided for boxplot visualization")
           do.call(boxplot, c(list(formula = as.formula(paste(y, "~", x)), data = data, main = main, xlab = xlab, ylab = ylab), params))
         }
  )
}
