#' Simulate simultaneous power for multiple tests in a single model
#'
#' Simulate power to simultaneously detect significant effects for a set of statistical tests in a single multivariate model. This function simulates data based on a correlation matrix imposed using the corrvar function from the SimCorrMix package (Fialkowski, 2018) using Fleishman's third-order polynomial transformation (Fleishman, 1978), and can be used to estimate power for multivariate models with between 2 and 10 predictor variables and a single dependent variable.
#'
#' When you use this function (and we hope you do!), please cite it as:
#'
#' Le Forestier, J. M., Page-Gould, E., & Chasteen, A. L. (Forthcoming). Statistical power for a set of tests.
#'
#' @usage simusim.multivar(n = NULL, es = NULL, es1 = NULL, es2 = NULL)
#'
#' @param n Set the size of each sample to be drawn from the population. This is the sample size for which you are estimating statistical power. In other words, setting n to equal 100 will estimate statistical power at n = 100. Accepts any positive number smaller than your population.
#' @param es Set the units in which you are specifying your effect sizes. Accepts "d" for Cohen's d, "r" for correlation coefficients, and "r2" for percent of variance accounted for.
#' @param predictors How many predictor variables would you like to simulate and ultimately calculate power for? Default = 2. Accepts whole numbers in the range of 2 to 10. Note that this argument is required if you specify effect size values for more that 2 predictors.
#' @param popsize What is the size of the population you would like to simulate? This is the population from which you will ultimately draw your samples. Default = 100,000. Accepts any positive whole number
#' @param iterations How many times you would like to estimate your model in random samples drawn from your population? One model will be run in each random sample. Default = 5,000. Accepts any whole number greater than 0.
#' @param alpha Set your alpha level. This is the threshold below which p-values will be considered significant. Default = 0.05. Accepts any number greater than 0 and less than 1.
#' @param seed Set a seed to make your results reproducible. Default = 1. Accepts any number.
#' @param es1...es10 The effect size, expressed in units specified in the es argument, for the relationship between each predictor and the dependent variable. You should specify a number of "es"s equal to the number you specified in the "predictors" argument. That is, if you set predictors to equal 4, you should supply values for es1, es2, es3, and es4. You should always specify these in order, beginning with es1, and not skipping any. Accepts any number.
#' @param iv1iv2_cov...iv9iv10_cov The relationships between each set of predictors, specified in the same units as your effect sizes. Specifying relationships between predictors is optional unless your predictors, together, account for more than 100% of the variance in your DV, in which case you must specify relationships between your predictors to make that possible. Note that not specifying relationships between predictors in the same model may lead to overestimation of statistical power, particularly in the presence of large effect sizes. Default = 0. Accepts any number.
#'
#' @return A dataframe containing a power estiamte, expressed as a decimal, for each of the effects individually, and for all the effects simultaneously.
#'
#' @author Joel Le Forestier (joel.leforestier@@mail.utoronto.ca)
#'
#' @references Fialkowski, A. C. (2018). SimmCorrMix: Simulation of correlated data with multiple variable types including continuous and count mixture distributions. Comprehensive R Archive Network (CRAN).
#'
#' Fleishman, A. I. (1978). A method for simulating non-normal distributions. Psychometrika, 43, 521-532.
#'
#' Le Forestier, J. M., Page-Gould, E., & Chasteen, A. (Forthcoming). Statistical power for a set of tests.
#'
#' @examples # A basic example, leaving all the defaults in place.
#'
#' simusim.multivar(n = 150, es = "r", es1 = .2, es2 = .15)
#'
#' # Another example, customizing additional parameters.
#'
#' simusim.multivar(n = 300, es = "d", es1 = .4, es2 = .3, es3 = .45, predictors = 3,
#'      iv1iv2_cov = .2, popsize = 500000, alpha = .01, seed = 123)
#'
#' @export

simusim.multivar <- function(n, es, es1, es2,
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
                    predictors = 2, popsize = 100000, iterations = 5000, alpha = .05, seed = 1) {

  # Throw a warning if the user has specified the wrong number of predictors #
  dummy_beta <- 0


  # Throw a warning if the user has specified a number of parameters not equal to the number of effect sizes they specified #
  specified_params <- table(c(es1, es2, es3, es4, es5, es6, es7, es8, es9, es10, dummy_beta))

  if(specified_params["0"] != (11 - predictors)) {
    stop("You have specified the wrong number of predictors. Make sure to set the predictors argument to equal the number of predictors for which you have specified effect sizes. Additionally, ensure to specify effect sizes for predictors in order, beginning with es1, and not skipping any along the way.")
  }

  # Throw a warning if the user has left n blank #
  if(missing(n)) {
    stop("Missing value for n. Please specify a sample size for your tests.")
  }

  # Throw a warning if the sample size is greater than or equal to the population size #
  if(n >= popsize) {
    stop("Your sample size is greater than or equal to your population size. Please increase popsize or decrease n.")
  }

  # Throw a warning if the number of iterations is 0 or negative #
  if(iterations < 1 | round(iterations, 0) != iterations) {
    stop("You have specified an invalid number of iterations. Please specify a whole number greater than 0.")
  }

  # Throw a warning if the sample size is 0 or negative #
  if(n < 1 | round(n, 0) != n) {
    stop("You have specified an invalid sample size. Please specify a whole number greater than 0.")
  }

  # Throw a warning if the population size is 0 or negative #
  if(popsize < 1 | round(popsize, 0) != popsize) {
    stop("You have specified an invalid population size. Please specify a whole number greater than 0.")
  }

  # Throw a warning if es is not r, r2, or d #
  if(es != "d" & es != "r" & es != "r2") {
    stop("You have specified an invalid type effect size. Please specify whether you are inputting")
  }

  # Throw a warning if the alpha level is not a number > 0 and < 1 #
  if(alpha <=0 | alpha >= 1) {
    stop("You have specified an incorrect alpha level. Please specify an alpha level greater than 0 and less than 1.")
  }

  # Generate the correlation matrix for the population #
  cortable <- matrix(c(1, iv1iv2_cov, iv1iv3_cov, iv1iv4_cov, iv1iv5_cov, iv1iv6_cov, iv1iv7_cov, iv1iv8_cov, iv1iv9_cov, iv1iv10_cov, es1,
                       iv1iv2_cov, 1, iv2iv3_cov, iv2iv4_cov, iv2iv5_cov, iv2iv6_cov, iv2iv7_cov, iv2iv8_cov, iv2iv9_cov, iv2iv10_cov, es2,
                       iv1iv3_cov, iv2iv3_cov, 1, iv3iv4_cov, iv3iv5_cov, iv3iv6_cov, iv3iv7_cov, iv3iv8_cov, iv3iv9_cov, iv3iv10_cov, es3,
                       iv1iv4_cov, iv2iv4_cov, iv3iv4_cov, 1, iv4iv5_cov, iv4iv6_cov, iv4iv7_cov, iv4iv8_cov, iv4iv9_cov, iv4iv10_cov, es4,
                       iv1iv5_cov, iv2iv5_cov, iv3iv5_cov, iv4iv5_cov, 1, iv5iv6_cov, iv5iv7_cov, iv5iv8_cov, iv5iv9_cov, iv5iv10_cov, es5,
                       iv1iv6_cov, iv2iv6_cov, iv3iv6_cov, iv4iv6_cov, iv5iv6_cov, 1, iv6iv7_cov, iv6iv8_cov, iv6iv9_cov, iv6iv10_cov, es6,
                       iv1iv7_cov, iv2iv7_cov, iv3iv7_cov, iv4iv7_cov, iv5iv7_cov, iv6iv7_cov, 1, iv7iv8_cov, iv7iv9_cov, iv7iv10_cov, es7,
                       iv1iv8_cov, iv2iv8_cov, iv3iv8_cov, iv4iv8_cov, iv5iv8_cov, iv6iv8_cov, iv7iv8_cov, 1, iv8iv9_cov, iv8iv10_cov, es8,
                       iv1iv9_cov, iv2iv9_cov, iv3iv9_cov, iv4iv9_cov, iv5iv9_cov, iv6iv9_cov, iv7iv9_cov, iv8iv9_cov, 1, iv9iv10_cov, es9,
                       iv1iv10_cov, iv2iv10_cov, iv3iv10_cov, iv4iv10_cov, iv5iv10_cov, iv6iv10_cov, iv7iv10_cov, iv8iv10_cov, iv9iv10_cov, 1, es10,
                       es1, es2, es3, es4, es5, es6, es7, es8, es9, es10, 1), 11, 11)

  # Convert effect size into correlation coefficients #
  if (es == "d") {
    cortable <- cortable / sqrt(cortable**2 + 4)
    for(a in 1:11) {
    cortable[a,a] <- 1
    }
  } else if (es == "r2") {
    cortable <- sqrt(cortable)
    for(a in 1:11) {
      cortable[a,a] <- 1
    }
  }

  # Let the user know it's working #
  message("Simulating the population.")

  # Simulate the population #
  set.seed(seed)
  invisible(capture.output(dist <- SimCorrMix::corrvar(n = popsize,
                   k_cont = 11,
                   method = "Fleishman",
                   means = rep(0, times = 11),
                   vars = rep(1, times = 11),
                   skews = rep(0, times = 11),
                   skurts = rep(0, times = 11),
                   rho = cortable)$Y_cont))

  dist <- data.frame(dist)

  names(dist) <- c("iv1", "iv2", "iv3", "iv4", "iv5", "iv6", "iv7", "iv8", "iv9", "iv10", "dv")

  # Set up the model #
  design <- "dv ~ iv1"
  for(v in 2:predictors) {
    design <- paste(design, " + iv", v, sep = "")
  }

  # Let the user know it's working #
  message(paste("Running", iterations, "sets of tests. This may take a minute.", sep = " "))

  # Set up for parallelization #
  doParallel::registerDoParallel(parallel::detectCores())
  `%dopar%` <- foreach::`%dopar%`

  # Run the model i times #
  result <- vector()

  result <- foreach::foreach (i=1:iterations, .combine=rbind) %dopar% {
    set.seed(seed + i)

    sample <- dist[sample(nrow(dist), size = n),]
    model <- data.frame(summary(lm(design, data = sample))$coefficients)

    for(v in 1:predictors){
      assign(paste("es", v, "_result", sep = ""), model$Pr...t..[v+1] < alpha)
    }

    if (predictors == 2) {
      result <- rbind(result, data.frame(es1_result, es2_result))
    } else if (predictors == 3) {
      result <- rbind(result, data.frame(es1_result, es2_result, es3_result))
    } else if (predictors == 4) {
      result <- rbind(result, data.frame(es1_result, es2_result, es3_result, es4_result))
    } else if (predictors == 5) {
      result <- rbind(result, data.frame(es1_result, es2_result, es3_result, es4_result, es5_result))
    } else if (predictors == 6) {
      result <- rbind(result, data.frame(es1_result, es2_result, es3_result, es4_result, es5_result, es6_result))
    } else if (predictors == 7) {
      result <- rbind(result, data.frame(es1_result, es2_result, es3_result, es4_result, es5_result, es6_result, es7_result))
    } else if (predictors == 8) {
      result <- rbind(result, data.frame(es1_result, es2_result, es3_result, es4_result, es5_result, es6_result, es7_result, es8_result))
    } else if (predictors == 9) {
      result <- rbind(result, data.frame(es1_result, es2_result, es3_result, es4_result, es5_result, es6_result, es7_result, es8_result, es9_result))
    } else if (predictors == 10) {
      result <- rbind(result, data.frame(es1_result, es2_result, es3_result, es4_result, es5_result, es6_result, es7_result, es8_result, es9_result, es10_result))
    }

  }

  # Create variable for calculating simultaneous power #
  result$simultaneous <- FALSE

  if (predictors == 2) {
    result$simultaneous[result$es1_result == TRUE & result$es2_result == TRUE] <- TRUE
  } else if (predictors == 3) {
    result$simultaneous[result$es1_result == TRUE & result$es2_result == TRUE & result$es3_result == TRUE] <- TRUE
  } else if (predictors == 4) {
    result$simultaneous[result$es1_result == TRUE & result$es2_result == TRUE & result$es3_result == TRUE & result$es4_result == TRUE] <- TRUE
  } else if (predictors == 5) {
    result$simultaneous[result$es1_result == TRUE & result$es2_result == TRUE & result$es3_result == TRUE & result$es4_result == TRUE & result$es5_result == TRUE] <- TRUE
  } else if (predictors == 6) {
    result$simultaneous[result$es1_result == TRUE & result$es2_result == TRUE & result$es3_result == TRUE & result$es4_result == TRUE & result$es5_result == TRUE &
                   result$es6_result == TRUE] <- TRUE
  } else if (predictors == 7) {
    result$simultaneous[result$es1_result == TRUE & result$es2_result == TRUE & result$es3_result == TRUE & result$es4_result == TRUE & result$es5_result == TRUE &
                   result$es6_result == TRUE & result$es7_result == TRUE] <- TRUE
  } else if (predictors == 8) {
    result$simultaneous[result$es1_result == TRUE & result$es2_result == TRUE & result$es3_result == TRUE & result$es4_result == TRUE & result$es5_result == TRUE &
                   result$es6_result == TRUE & result$es7_result == TRUE & result$es8_result == TRUE] <- TRUE
  } else if (predictors == 9) {
    result$simultaneous[result$es1_result == TRUE & result$es2_result == TRUE & result$es3_result == TRUE & result$es4_result == TRUE & result$es5_result == TRUE &
                   result$es6_result == TRUE & result$es7_result == TRUE & result$es8_result == TRUE & result$es9_result == TRUE] <- TRUE
  } else if (predictors == 10) {
    result$simultaneous[result$es1_result == TRUE & result$es2_result == TRUE & result$es3_result == TRUE & result$es4_result == TRUE & result$es5_result == TRUE &
                   result$es6_result == TRUE & result$es7_result == TRUE & result$es8_result == TRUE & result$es9_result == TRUE & result$es10_result == TRUE] <- TRUE
  }

  result[result == "TRUE"] <- 1

  # Calculate power #
  simultaneous_power <- mean(result$simultaneous)
  es1_power <- mean(result$es1_result)
  es2_power <- mean(result$es2_result)
  if (predictors > 2) {
    es3_power <- mean(result$es3_result)
  }
  if (predictors > 3) {
    es4_power <- mean(result$es4_result)
  }
  if (predictors > 4) {
    es5_power <- mean(result$es5_result)
  }
  if (predictors > 5) {
    es6_power <- mean(result$es6_result)
  }
  if (predictors > 6) {
    es7_power <- mean(result$es7_result)
  }
  if (predictors > 7) {
    es8_power <- mean(result$es8_result)
  }
  if (predictors > 8) {
    es9_power <- mean(result$es9_result)
  }
  if (predictors > 9) {
    es10_power <- mean(result$es10_result)
  }

  # Summarize results

  if (predictors == 2) {
    power <- data.frame(es1_power, es2_power, simultaneous_power)
  } else if ( predictors == 3) {
    power <- data.frame(es1_power, es2_power, es3_power, simultaneous_power)
  } else if ( predictors == 4) {
    power <- data.frame(es1_power, es2_power, es3_power, es4_power, simultaneous_power)
  } else if ( predictors == 5) {
    power <- data.frame(es1_power, es2_power, es3_power, es4_power, es5_power, simultaneous_power)
  } else if ( predictors == 6) {
    power <- data.frame(es1_power, es2_power, es3_power, es4_power, es5_power, es6_power, simultaneous_power)
  } else if ( predictors == 7) {
    power <- data.frame(es1_power, es2_power, es3_power, es4_power, es5_power, es6_power, es7_power, simultaneous_power)
  } else if ( predictors == 8) {
    power <- data.frame(es1_power, es2_power, es3_power, es4_power, es5_power, es6_power, es7_power, es8_power, simultaneous_power)
  } else if ( predictors == 9) {
    power <- data.frame(es1_power, es2_power, es3_power, es4_power, es5_power, es6_power, es7_power, es8_power, es9_power, simultaneous_power)
  } else if ( predictors == 10) {
    power <- data.frame(es1_power, es2_power, es3_power, es4_power, es5_power, es6_power, es7_power, es8_power, es9_power, es10_power, simultaneous_power)
  }

  return(power)
}
