#' Simulate and visualize simultaneous power curves
#'
#' Simulate and visualize individual and simultaneous power at a range of sample sizes of a set of statistical tests. Simultaneous power simulations are conducted using \link[SimuSim]{simusim.multivars} or \link[SimuSim]{simusim.multimodels}. A detailed walkthrough and set of vignettes for this and other SimuSim functions is available [here.](https://doi.org/10.31219/osf.io/w96uk)
#'
#' The power.curves function is the first step in the suggested SimuSim workflow. Start here to visualize the approximate simultaneous power space occupied by your set of tests, then use either \link[SimuSim]{simusim.multivars} or \link[SimuSim]{simusim.multimodels} with lrager numbers of iterations for final power calculations with more stable estimates.
#'
#' When you use this function (and we hope you do!), please cite the package:
#'
#' Le Forestier, J. M. (2020). SimuSim: Simultaneous power analysis for a set of statistical tests. https://doi.org/10.31219/osf.io/w96uk
#'
#' and/or cite the accompanying paper:
#'
#' Le Forestier, J. M., Page-Gould, E., & Chasteen, A. L. (Forthcoming). Statistical power for a set of tests.
#'
#' @usage power.curve(method = NULL, min = NULL, max = NULL, increment = 20, thresholds = c(80, 90, 95), es_units = NULL, es1 = NULL, es2 = NULL, es3...es10 = 0, iv1iv2_cov...iv9iv10_cov = 0, predictors = 2, models = 2, null_effect = 0, popsize = 100000, iterations = 1000, alpha = .05, bonferroni = FALSE, seed = 1, iv1iv2_cov...iv9iv10_cov = 0)
#'
#' @param method Specify which SimuSim function you would like to use to computer simultaneous power. Accepts either "simusim.multivars" or "simusim.multimodels". This argument has no default.
#' @param min Set the minimum sample size for which you would like to estimate simultaneous power. Accepts any positive number smaller than your population and your maximum sample size. This argument has no default.
#' @param max Set the maximum sample size for which you would like to estimate simultaneous power. Accepts any positive number smaller than your population and larger than your maximum sample size. This argument has no default.
#' @param increment Set the increments by which you would like the sample size to increase between power analyses. Accepts any positive whole number. Default = 20.
#' @param thresholds Set the power thresholds for which you would like sample size estimates printed. Accepts either any single number between 0 and 100 or a vector of numbers between 0 and 100. Default = c(80, 90, 95).
#' @param es_units Set the units in which you are specifying your effect sizes. Accepts "d" for Cohen's d, "r" for correlation coefficients, and "r2" for percent of variance accounted for. This argument has no default.
#' @param predictors If you are estimating power for multiple tests in a single model, how many predictor variables would you like to simulate and ultimately calculate power for? Default = 2. Note that this argument is required if you specify effect size values for more that 2 predictors. Accepts whole numbers in the range of 2 to 10. Default = 2.
#' @param models If you are estimating power for tests in separate models, how many models would you like to simulate and ultimately calculate power for? Default = 2. Accepts whole numbers in the range of 2 to 10. Note that this argument is required if you specify effect size values for more that 2 predictors.
#' @param null_effect For which, if any, of your predictors or models are you computing "null power?" If you want to compute "power" to NOT detect an effect, use this argument to specify which effects are predicted nulls by setting this argument equal to the number(s) corresponding to the predictors you hypothesize to be null. Accepts either a single whole number between 1 and the number of predictors you have specified or a vector of numbers between 1 and the the number of predictors you have specified. Default = no null effects.
#' @param popsize What is the size of the population you would like to simulate? This is the population from which you will ultimately draw your samples. Note that the population you simulate does NOT have to be the same size as the real-world population to which you intend to generalize your results, and that simulating very large populations may require more computer memory than is available to some users. Accepts any positive whole number. Default = 100,000.
#' @param iterations How many times you would like to run your model in random samples drawn from your population? One model will be run in each random sample. Accepts any whole number greater than 0. Default = 1,000.
#' @param alpha Set your alpha level. This is the threshold below which p-values will be considered significant. Accepts any number greater than 0 and less than 1. Default = 0.05.
#' @param bonferroni Apply a bonferroni correction? This is suggested if you intend on interpreting the results of multiple tests individually, but not if you intend on assessing a single research question by triangulating across multiple tests (Le Forestier, Page-Gould, & Chasteen, Forthcoming). Accepts TRUE or FALSE. Default = FALSE.
#' @param seed Set a seed to make your results reproducible. Accepts any number. Default = 1.
#' @param es1...es10 The effect size, expressed in units specified in the es_units argument, for the relationship between each predictor and the dependent variable. You should specify a number of "es"s equal to the number you specified in the "predictors" argument. That is, if you set predictors to equal 4, you should supply values for es1, es2, es3, and es4. You should always specify these in order, beginning with es1, and not skipping any. Accepts any number. These arguments have no defaults.
#' @param iv1iv2_cov...iv9iv10_cov The relationships between each set of predictors, specified in correlation coefficients. Specifying relationships between predictors is optional unless your predictors, together, account for more than 100% of the variance in your DV, in which case you must specify relationships between your predictors to make that possible. Accepts any number between -1 and 1. Default = 0.
#'
#' @return A plot of power curves and a set of sample size estimates for default or user-defined power thresholds.
#'
#' A dataframe containing a power estiamte, expressed as a decimal, for each of the effects individually, and for all the effects simultaneously.
#'
#' @author Joel Le Forestier (joel.leforestier@@mail.utoronto.ca)
#'
#' @references Le Forestier, J. M., Page-Gould, E., & Chasteen, A. (Forthcoming). Statistical power for a set of tests.
#'
#' @examples # An example using simusim.multivars
#'
#' power.curve(method = "simusim.multivars", min = 200, max = 500, predictors = 3,
#' es_units = "d", es1 = .35, es2 = .45, es3 = .50)
#'
#' # An example using simusim.multimodels
#'
#' power.curve(method = "simusim.multimodels", min = 100, max = 400, increment = 40, models = 7,
#' es_units = "r2", es1 = .05, es2 = .06, es3 = .08, es4 = .06, es5 = .08, es6 = .07, es7 = .04)
#'
#' @export

power.curve <- function(method, min, max, increment = 20, thresholds = c(80, 90, 95),
                              es_units, es1, es2,
                              es3 = 0, es4 = 0, es5 = 0, es6 = 0, es7 = 0, es8 = 0, es9 = 0, es10 = 0,
                              iv1iv2_cov = 0, iv1iv3_cov = 0, iv1iv4_cov = 0, iv1iv5_cov = 0, iv1iv6_cov = 0, iv1iv7_cov = 0, iv1iv8_cov = 0, iv1iv9_cov = 0, iv1iv10_cov = 0,
                              iv2iv3_cov = 0, iv2iv4_cov = 0, iv2iv5_cov = 0, iv2iv6_cov = 0, iv2iv7_cov = 0, iv2iv8_cov = 0, iv2iv9_cov = 0, iv2iv10_cov = 0,
                              iv3iv4_cov = 0, iv3iv5_cov = 0, iv3iv6_cov = 0, iv3iv7_cov = 0, iv3iv8_cov = 0, iv3iv9_cov = 0, iv3iv10_cov = 0,
                              iv4iv5_cov = 0, iv4iv6_cov = 0, iv4iv7_cov = 0, iv4iv8_cov = 0, iv4iv9_cov = 0, iv4iv10_cov = 0,
                              iv5iv6_cov = 0, iv5iv7_cov = 0, iv5iv8_cov = 0, iv5iv9_cov = 0, iv5iv10_cov = 0,
                              iv6iv7_cov = 0, iv6iv8_cov = 0, iv6iv9_cov = 0, iv6iv10_cov = 0,
                              iv7iv8_cov = 0, iv7iv9_cov = 0, iv7iv10_cov = 0,
                              iv8iv9_cov = 0, iv8iv10_cov = 0,
                              iv9iv10_cov = 0,
                              predictors = 2, models = 2, null_effect = 0, popsize = 100000, iterations = 1000, alpha = .05, bonferroni = FALSE, seed = 1) {

  # Throw a warning if method is wrong or missing #
  if(method != "simusim.multivars" & method != "simusim.multimodels") {
    stop("You have not correctly specified which SimuSim function you would like to use to estimate simutaneous power.")
  }

  # Throw a warning if min is greater than or equal to the population size or max size#
  if(min >= popsize | min >= max) {
    stop("Your minimum sample size is greater than or equal to either your population size or your maximum sample size. Please adjust popsize, min, or max.")
  }

  # Throw a warning if max is greater than or equal to the population size or less than or equal to min size#
  if(max >= popsize | max <= min) {
    stop("Your maximum sample size is greater than or equal to either your population size or less than or equal to your minimum sample size. Please adjust popsize, min, or max.")
  }

  # Throw a warning if the increment is 0, negative, or not a whole number #
  if(increment < 1 | round(increment, 0) != increment) {
    stop("You have specified an invalid increment size. Please specify a whole number greater than 0.")
  }

  # Throw a warning if min or max sample size is 0, negative, or not a whole number #
  if(min < 1 | round(min, 0) != min | max < 1 | round(max, 0) != max) {
    stop("You have specified an invalid sample size. Please specify a whole number greater than 0.")
  }

  # Throw a warning if any of the thresholds are > 100 or < 0 #
  error_high <- thresholds > 100
  error_low <- thresholds < 0

  if(TRUE %in% error_high | TRUE %in% error_low) {
    stop("You have specified an impossible power threshold in the \"thresholds\" aregument. Please specify a number, or vector of numbers, between 0 and 1.")
  }

  # Select a function and perform the power analyses #
  if(method == "simusim.multivars") {

    data <- data.frame()

    # let them know it's working
    message(paste("Simulating", (max - min) / increment + 1, "populations and performing", iterations, "sets of tests in each. This may take a few minutes."))

    # run the power analyses #
    for (n in seq(from=min, to=max, by=increment)) {
      data <- rbind(data,
                    suppressMessages(SimuSim::simusim.multivars(
                      n = n,
                      es_units = es_units,
                      predictors = predictors,
                      null_effect = null_effect,
                      popsize = popsize,
                      iterations = iterations,
                      alpha = alpha,
                      seed = seed,
                      es1 = es1, es2 = es2, es3 = es3, es4 = es4, es5 = es5, es6 = es6, es7 = es7, es8 = es8, es9 = es9, es10 = es10,
                      iv1iv2_cov = iv1iv2_cov, iv1iv3_cov = iv1iv3_cov, iv1iv4_cov = iv1iv4_cov, iv1iv5_cov = iv1iv5_cov, iv1iv6_cov = iv1iv6_cov, iv1iv7_cov = iv1iv7_cov, iv1iv8_cov = iv1iv8_cov, iv1iv9_cov = iv1iv9_cov, iv1iv10_cov = iv1iv10_cov,
                      iv2iv3_cov = iv2iv3_cov, iv2iv4_cov = iv2iv4_cov, iv2iv5_cov = iv2iv5_cov, iv2iv6_cov = iv2iv6_cov, iv2iv7_cov = iv2iv7_cov, iv2iv8_cov = iv2iv8_cov, iv2iv9_cov = iv2iv9_cov, iv2iv10_cov = iv2iv10_cov,
                      iv3iv4_cov = iv3iv4_cov, iv3iv5_cov = iv3iv5_cov, iv3iv6_cov = iv3iv6_cov, iv3iv7_cov = iv3iv7_cov, iv3iv8_cov = iv3iv8_cov, iv3iv9_cov = iv3iv9_cov, iv3iv10_cov = iv3iv10_cov,
                      iv4iv5_cov = iv4iv5_cov, iv4iv6_cov = iv4iv6_cov, iv4iv7_cov = iv4iv7_cov, iv4iv8_cov = iv4iv8_cov, iv4iv9_cov = iv4iv9_cov, iv4iv10_cov = iv4iv10_cov,
                      iv5iv6_cov = iv5iv6_cov, iv5iv7_cov = iv5iv7_cov, iv5iv8_cov = iv5iv8_cov, iv5iv9_cov = iv5iv9_cov, iv5iv10_cov = iv5iv10_cov,
                      iv6iv7_cov = iv6iv7_cov, iv6iv8_cov = iv6iv8_cov, iv6iv9_cov = iv6iv9_cov, iv6iv10_cov = iv6iv10_cov,
                      iv7iv8_cov = iv7iv8_cov, iv7iv9_cov = iv7iv9_cov, iv7iv10_cov = iv7iv10_cov,
                      iv8iv9_cov = iv8iv9_cov, iv8iv10_cov = iv8iv10_cov,
                      iv9iv10_cov = iv9iv10_cov,
                      print_result = FALSE)))
    }

    data$n <- seq(from = min, to = max, by = increment)

    # plot the result #
    color_options <- c("#9F018A", #color blind-friendly palette
                       "#009F81",
                       "#FF5AAF",
                       "#00FCCF",
                       "#8400CD",
                       "#008DF9",
                       "#A40122",
                       "#E20134",
                       "#FF6E3A",
                       "#FFC33B")
    palette <- append(color_options[1:predictors], "#000000") #simultaneous power is always black

    longdata <- tidyr::pivot_longer(data = data,
                                    cols = -n,
                                    names_to = "param",
                                    values_to = "estimate")
    names(longdata) <- c("n", "parameter", "power")

    for(t in 1:predictors) {
      longdata$parameter[longdata$parameter == paste0("es", t, "_power")] <- paste("predictor", t)
    }
    longdata$parameter[longdata$parameter == "simultaneous_power"] <- "simultaneous power"

    plot <- ggplot2::ggplot(data = longdata, mapping = ggplot2::aes(x = n, y = power)) +
      ggplot2::scale_color_manual(name = "Parameter",
                                  values = palette) +
      ggplot2::geom_line(ggplot2::aes(color = parameter), size = 1) +
      ggplot2::scale_y_continuous(limits = c(0, 100),
                         breaks=c(0, 20, 40, 60, 80, 100),
                         labels=c("0%", "20%", "40%", "60%", "80%", "100%")) +
      ggplot2::xlab("Sample Size") +
      ggplot2::ylab("Power") +
      ggplot2::theme_minimal()

  } else if (method == "simusim.multimodels") {

    data <- data.frame()

    # let them know it's working
    message(paste("Simulating", (max - min) / increment + 1, "populations and performing", iterations, "sets of tests in each. This may take a few minutes."))

    # run the power analyses #
    for (n in seq(from=min, to=max, by=increment)) {
      data <- rbind(data,
                    suppressMessages(SimuSim::simusim.multimodels(
                      n = n,
                      es_units = es_units,
                      models = models,
                      null_effect = null_effect,
                      popsize = popsize,
                      iterations = iterations,
                      alpha = alpha,
                      seed = seed,
                      es1 = es1, es2 = es2, es3 = es3, es4 = es4, es5 = es5, es6 = es6, es7 = es7, es8 = es8, es9 = es9, es10 = es10,
                      print_result = FALSE)))
    }

  data$n <- seq(from = min, to = max, by = increment)

    # plot the result #
    color_options <- c("#9F018A", #color blind-friendly palette
                       "#009F81",
                       "#FF5AAF",
                       "#00FCCF",
                       "#8400CD",
                       "#008DF9",
                       "#A40122",
                       "#E20134",
                       "#FF6E3A",
                       "#FFC33B")
    palette <- append(color_options[1:models], "#000000") #simultaneous power is always black

    longdata <- tidyr::pivot_longer(data = data,
                                    cols = -n,
                                    names_to = "param",
                                    values_to = "estimate")
    names(longdata) <- c("n", "model", "power")

    for(t in 1:models) {
      longdata$model[longdata$model == paste0("es", t, "_power")] <- paste("model", t)
    }
    longdata$model[longdata$model == "simultaneous_power"] <- "simultaneous power"

    plot <- ggplot2::ggplot(data = longdata, mapping = ggplot2::aes(x = n, y = power)) +
      ggplot2::scale_color_manual(name = "Model",
                                  values = palette) +
      ggplot2::geom_line(ggplot2::aes(color = model), size = 1) +
      ggplot2::scale_y_continuous(limits = c(0, 100),
                                  breaks=c(0, 20, 40, 60, 80, 100),
                                  labels=c("0%", "20%", "40%", "60%", "80%", "100%")) +
      ggplot2::xlab("Sample Size") +
      ggplot2::ylab("Power") +
      ggplot2::theme_minimal()
  }

  # print the sample size estimates #
  thresh_n <- list()
  for (x in 1:length(thresholds)) {
    if (thresholds[x] < min(data$simultaneous_power)) {
      thresh_n <- append(thresh_n, NA)
      thresh_n <- append(thresh_n, NA)
    } else {
      thresh_n <- append(thresh_n, data$n[which(data$simultaneous_power >= thresholds[x])[1] - 1])
      thresh_n <- append(thresh_n, data$n[which(data$simultaneous_power >= thresholds[x])[1]])
    }
    if (is.na(thresh_n[[x * 2 - 1]]) | is.na(thresh_n[[x * 2]])) {
      message(paste0("The ", thresholds[x], "% simultaneous power threshold was not crossed."))
    } else {
      cat("The ", thresholds[x], "% simultaneous power threshold was crossed between n = ", thresh_n[[x * 2 - 1]], " and n = ", thresh_n[[x * 2]], "\n", sep = "")
    }
  }

  # return the plot #
  return(plot)
}



