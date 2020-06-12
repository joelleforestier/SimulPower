#' Simulate simultaneous power for multiple tests
#'
#' Simulate power to simultaneously detect significant effects for a set of statistical tests in a single multivariate model. This function simulates data based on a correlation matrix imposed using the SimMultiCorrData package (Fialkowski, A. C., 2018) using Fleishman's third-order polynomial transformation (Fleishman, 1978), and can be used to estimate power for multivariate models with between 2 and 10 predictor variables and a single dependent variable.
#'
#' @usage simusim(n = NULL, b1 = NULL, b2 = NULL)
#'
#' @param predictors How many predictor variables would you like to simulate and ultimately calculate power for? Default = 2. Accepts integers in the range of 2 to 10.
#' @param popsize What is the size of the population you would like to simulate? This is the population from which you will ultimately draw you samples. Default = 100,000. Accepts any positive integer.
#' @param iterations How many times you would like to estimate your models in random samples drawn from your population? One model will be run in each random sample. Default = 10,000. Accepts any positive integer.
#' @param alpha Set your alpha level. This is the threshold below which p-values will be considered significant. Default = 0.05. Accepts any number.
#' @param seed Set a seed to make your results reproducible. Default = 1. Accepts any number.
#' @param n Set the size of each sample to be drawn from the population. This is the sample size for which you are estimating statistical power. In other words, setting n to equal 100 will estimate statistical power at n = 100.
#' @param b1...b10 The effect size, expressed as correlation coefficients (r), for each predictor in your model. You should specify a number of "b"s equal to the number you specified in the "predictors" argument. That is, if you set predictors to equal 4, you should supple values for b1, b2, b3, and b4. You should always specify these in order, beginning with b1, and not skipping any. You must specify effect sizes for at least b1 and b2, neither of which has a default value. Default values for all other bs are 0. Leave these defaults in place for any unused bs. Accepts any number between -1 and 1.
#' @param iv1iv2_cov...iv9iv10_cov The covariance, expressed as correlation coefficients (r), between each set of predictors. Specifying covariances between predictors is optional unless your predictors, together, account for more than 100% of the variance in your DV, in which case you must specify covariances between your predictors to make that possible. Default = 0. Accepts any number between -1 and 1.
#'
#' @value A dataframe containing a power estiamte, expressed as a decimal, for each of the effects individually, and for all the effects simultaneously (labelled "total_power").
#'
#' @references Fialkowski, A. C. (2018). SimMultiCorrData: Simulation of correlated data with multiple variable types. Comprehensive R Archive Network (CRAN).
#'
#' Fleishman, A. I. (1978). A method for simulating non-normal distributions. Psychometrika, 43, 521-532.
#'
#' Le Forestier, J. M., Page-Gould, E., & Chasteen, A. (Forthcoming). Statistical power for a set of tests.
#'
#' @examples
#' # A basic example, leaving all the defaults in place.
#'
#' simusim(n = 150, b1 = .2, b2 = .15)
#'
#' # Another example, customizing additional parameters.
#'
#' simusim(n = 300, b1 = .7, b2 = .6, b3 = .5, iv1iv2_cov = .5,
#'      popsize = 1000000, alpha = .01, seed = 12345)
#'
#' @export

simusim <- function(predictors = 2, popsize = 100000, iterations = 10000, alpha = .05, seed = 1,
                    n, b1, b2,
                    b3 = 0, b4 = 0, b5 = 0, b6 = 0, b7 = 0, b8 = 0, b9 = 0, b10 = 0,
                    iv1iv2_cov = 0, iv1iv3_cov = 0, iv1iv4_cov = 0, iv1iv5_cov = 0, iv1iv6_cov = 0, iv1iv7_cov = 0, iv1iv8_cov = 0, iv1iv9_cov = 0, iv1iv10_cov = 0,
                    iv2iv3_cov = 0, iv2iv4_cov = 0, iv2iv5_cov = 0, iv2iv6_cov = 0, iv2iv7_cov = 0, iv2iv8_cov = 0, iv2iv9_cov = 0, iv2iv10_cov = 0,
                    iv3iv4_cov = 0, iv3iv5_cov = 0, iv3iv6_cov = 0, iv3iv7_cov = 0, iv3iv8_cov = 0, iv3iv9_cov = 0, iv3iv10_cov = 0,
                    iv4iv5_cov = 0, iv4iv6_cov = 0, iv4iv7_cov = 0, iv4iv8_cov = 0, iv4iv9_cov = 0, iv4iv10_cov = 0,
                    iv5iv6_cov = 0, iv5iv7_cov = 0, iv5iv8_cov = 0, iv5iv9_cov = 0, iv5iv10_cov = 0,
                    iv6iv7_cov = 0, iv6iv8_cov = 0, iv6iv9_cov = 0, iv6iv10_cov = 0,
                    iv7iv8_cov = 0, iv7iv9_cov = 0, iv7iv10_cov = 0,
                    iv8iv9_cov = 0, iv8iv10_cov = 0,
                    iv9iv10_cov = 0) {

  # Throw a warning if the user has specified the wrong number of predictors #
  dummy_beta <- 0

  specified_params <- table(c(b1, b2, b3, b4, b5, b6, b7, b8, b9, b10, dummy_beta))

  if(specified_params["0"] != (11 - predictors)) {
    stop("You have specified the wrong number of predictors. Make sure to set the predictors argument to equal the number of predictors for which you have specified effect sizes. Additionally, ensure to specify effect sizes for preductors in order, beginning with b1, and not skipping any along the way. ")
  }

  # Throw a warning if the user has left n blank #
  if(missing(n)) {
    stop("Missing value for n. Please specify a sample size for your tests.")
  }

  # Throw a warning if the sample size is greater than or equal to the population size #
  if(n >= popsize) {
    stop("Your sample size is greater than or equal to your population size. Please increase popsize or decrease n.")
  }

  # Simulate the population #
  cortable <- matrix(c(1, iv1iv2_cov, iv1iv3_cov, iv1iv4_cov, iv1iv5_cov, iv1iv6_cov, iv1iv7_cov, iv1iv8_cov, iv1iv9_cov, iv1iv10_cov, b1,
                       iv1iv2_cov, 1, iv2iv3_cov, iv2iv4_cov, iv2iv5_cov, iv2iv6_cov, iv2iv7_cov, iv2iv8_cov, iv2iv9_cov, iv2iv10_cov, b2,
                       iv1iv3_cov, iv2iv3_cov, 1, iv3iv4_cov, iv3iv5_cov, iv3iv6_cov, iv3iv7_cov, iv3iv8_cov, iv3iv9_cov, iv3iv10_cov, b3,
                       iv1iv4_cov, iv2iv4_cov, iv3iv4_cov, 1, iv4iv5_cov, iv4iv6_cov, iv4iv7_cov, iv4iv8_cov, iv4iv9_cov, iv4iv10_cov, b4,
                       iv1iv5_cov, iv2iv5_cov, iv3iv5_cov, iv4iv5_cov, 1, iv5iv6_cov, iv5iv7_cov, iv5iv8_cov, iv5iv9_cov, iv5iv10_cov, b5,
                       iv1iv6_cov, iv2iv6_cov, iv3iv6_cov, iv4iv6_cov, iv5iv6_cov, 1, iv6iv7_cov, iv6iv8_cov, iv6iv9_cov, iv6iv10_cov, b6,
                       iv1iv7_cov, iv2iv7_cov, iv3iv7_cov, iv4iv7_cov, iv5iv7_cov, iv6iv7_cov, 1, iv7iv8_cov, iv7iv9_cov, iv7iv10_cov, b7,
                       iv1iv8_cov, iv2iv8_cov, iv3iv8_cov, iv4iv8_cov, iv5iv8_cov, iv6iv8_cov, iv7iv8_cov, 1, iv8iv9_cov, iv8iv10_cov, b8,
                       iv1iv9_cov, iv2iv9_cov, iv3iv9_cov, iv4iv9_cov, iv5iv9_cov, iv6iv9_cov, iv7iv9_cov, iv8iv9_cov, 1, iv9iv10_cov, b9,
                       iv1iv10_cov, iv2iv10_cov, iv3iv10_cov, iv4iv10_cov, iv5iv10_cov, iv6iv10_cov, iv7iv10_cov, iv8iv10_cov, iv9iv10_cov, 1, b10,
                       b1, b2, b3, b4, b5, b6, b7, b8, b9, b10, 1), 11, 11)

  sink(file = "sink.txt")
  set.seed(seed)
  dist <- SimMultiCorrData::rcorrvar(n = popsize,
                   k_cont = 11,
                   method = "Fleishman",
                   means = rep(0, times = 11),
                   vars = rep(1, times = 11),
                   skews = rep(0, times = 11),
                   skurts = rep(0, times = 11),
                   rho = cortable)
  dist <- dist$continuous_variables
  sink(file = NULL)
  file.remove("sink.txt")

  names(dist) <- c("iv1", "iv2", "iv3", "iv4", "iv5", "iv6", "iv7", "iv8", "iv9", "iv10", "dv")

  # Set up the model #
  design <- "dv ~ iv1"
  for(v in 2:predictors) {
    design <- paste(design, " + iv", v, sep = "")
  }

  # Let the user know it's working #
  message(paste("Running", iterations, "simulations. This may take a minute.", sep = " "))

  # Run the model i times #
  result <- vector()

  set.seed(seed + 1)

  for(i in 1:iterations) {
    sample <- dist[sample(nrow(dist), size = n),]

    model <- data.frame(summary(lm(design, data = sample))$coefficients)

    for(v in 1:predictors){
      assign(paste("b", v, "_result", sep = ""), model$Pr...t..[v+1] < alpha)
    }

    if (predictors == 2) {
      result <- rbind(result, data.frame(b1_result, b2_result))
    } else if (predictors == 3) {
      result <- rbind(result, data.frame(b1_result, b2_result, b3_result))
    } else if (predictors == 4) {
      result <- rbind(result, data.frame(b1_result, b2_result, b3_result, b4_result))
    } else if (predictors == 5) {
      result <- rbind(result, data.frame(b1_result, b2_result, b3_result, b4_result, b5_result))
    } else if (predictors == 6) {
      result <- rbind(result, data.frame(b1_result, b2_result, b3_result, b4_result, b5_result, b6_result))
    } else if (predictors == 7) {
      result <- rbind(result, data.frame(b1_result, b2_result, b3_result, b4_result, b5_result, b6_result, b7_result))
    } else if (predictors == 8) {
      result <- rbind(result, data.frame(b1_result, b2_result, b3_result, b4_result, b5_result, b6_result, b7_result, b8_result))
    } else if (predictors == 9) {
      result <- rbind(result, data.frame(b1_result, b2_result, b3_result, b4_result, b5_result, b6_result, b7_result, b8_result, b9_result))
    } else if (predictors == 10) {
      result <- rbind(result, data.frame(b1_result, b2_result, b3_result, b4_result, b5_result, b6_result, b7_result, b8_result, b9_result, b10_result))
    }
  }

  # Create variable for calculating simultaneous power #
  result$total <- FALSE

  if (predictors == 2) {
    result$total[result$b1_result == TRUE & result$b2_result == TRUE] <- TRUE
  } else if (predictors == 3) {
    result$total[result$b1_result == TRUE & result$b2_result == TRUE & result$b3_result == TRUE] <- TRUE
  } else if (predictors == 4) {
    result$total[result$b1_result == TRUE & result$b2_result == TRUE & result$b3_result == TRUE & result$b4_result == TRUE] <- TRUE
  } else if (predictors == 5) {
    result$total[result$b1_result == TRUE & result$b2_result == TRUE & result$b3_result == TRUE & result$b4_result == TRUE & result$b5_result == TRUE] <- TRUE
  } else if (predictors == 6) {
    result$total[result$b1_result == TRUE & result$b2_result == TRUE & result$b3_result == TRUE & result$b4_result == TRUE & result$b5_result == TRUE &
                   result$b6_result == TRUE] <- TRUE
  } else if (predictors == 7) {
    result$total[result$b1_result == TRUE & result$b2_result == TRUE & result$b3_result == TRUE & result$b4_result == TRUE & result$b5_result == TRUE &
                   result$b6_result == TRUE & result$b7_result == TRUE] <- TRUE
  } else if (predictors == 8) {
    result$total[result$b1_result == TRUE & result$b2_result == TRUE & result$b3_result == TRUE & result$b4_result == TRUE & result$b5_result == TRUE &
                   result$b6_result == TRUE & result$b7_result == TRUE & result$b8_result == TRUE] <- TRUE
  } else if (predictors == 9) {
    result$total[result$b1_result == TRUE & result$b2_result == TRUE & result$b3_result == TRUE & result$b4_result == TRUE & result$b5_result == TRUE &
                   result$b6_result == TRUE & result$b7_result == TRUE & result$b8_result == TRUE & result$b9_result == TRUE] <- TRUE
  } else if (predictors == 10) {
    result$total[result$b1_result == TRUE & result$b2_result == TRUE & result$b3_result == TRUE & result$b4_result == TRUE & result$b5_result == TRUE &
                   result$b6_result == TRUE & result$b7_result == TRUE & result$b8_result == TRUE & result$b9_result == TRUE & result$b10_result == TRUE] <- TRUE
  }

  result[result == "TRUE"] <- 1

  # Calculate power #
  total_power <- mean(result$total)
  b1_power <- mean(result$b1_result)
  b2_power <- mean(result$b2_result)
  if (predictors > 2) {
    b3_power <- mean(result$b3_result)
  }
  if (predictors > 3) {
    b4_power <- mean(result$b4_result)
  }
  if (predictors > 4) {
    b5_power <- mean(result$b5_result)
  }
  if (predictors > 5) {
    b6_power <- mean(result$b6_result)
  }
  if (predictors > 6) {
    b7_power <- mean(result$b7_result)
  }
  if (predictors > 7) {
    b8_power <- mean(result$b8_result)
  }
  if (predictors > 8) {
    b9_power <- mean(result$b9_result)
  }
  if (predictors > 9) {
    b10_power <- mean(result$b10_result)
  }

  # Summarize results

  if (predictors == 2) {
    power <- data.frame(b1_power, b2_power, total_power)
  } else if ( predictors == 3) {
    power <- data.frame(b1_power, b2_power, b3_power, total_power)
  } else if ( predictors == 4) {
    power <- data.frame(b1_power, b2_power, b3_power, b4_power, total_power)
  } else if ( predictors == 5) {
    power <- data.frame(b1_power, b2_power, b3_power, b4_power, b5_power, total_power)
  } else if ( predictors == 6) {
    power <- data.frame(b1_power, b2_power, b3_power, b4_power, b5_power, b6_power, total_power)
  } else if ( predictors == 7) {
    power <- data.frame(b1_power, b2_power, b3_power, b4_power, b5_power, b6_power, b7_power, total_power)
  } else if ( predictors == 8) {
    power <- data.frame(b1_power, b2_power, b3_power, b4_power, b5_power, b6_power, b7_power, b8_power, total_power)
  } else if ( predictors == 9) {
    power <- data.frame(b1_power, b2_power, b3_power, b4_power, b5_power, b6_power, b7_power, b8_power, b9_power, total_power)
  } else if ( predictors == 10) {
    power <- data.frame(b1_power, b2_power, b3_power, b4_power, b5_power, b6_power, b7_power, b8_power, b9_power, b10_power, total_power)
  }

  return(power)
}
