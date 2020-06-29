#' Simulate simultaneous power for multiple tests in separate models
#'
#' Simulate power to simultaneously detect significant effects for a set of statistical tests in separate bivariate models. This function simulates data based on a correlation matrix imposed using the corrvar function from the SimCorrMix package (Fialkowski, 2018) using Fleishman's third-order polynomial transformation (Fleishman, 1978), and can be used to estimate power for up to 10 bivariate models with single independent variables and single dependent variables.
#'
#' When you use this function (and we hope you do!), please cite it as:
#'
#' Le Forestier, J. M., Page-Gould, E., & Chasteen, A. L. (Forthcoming). Statistical power for a set of tests.
#'
#' @usage simusim.bivar(n = NULL, b1 = NULL, b2 = NULL)
#'
#' @param models How many models would you like to simulate and ultimately calculate power for? Default = 2. Accepts whole numbers in the range of 2 to 10. Note that this argument is required if you specify more than 2 effect sizes. 
#' @param popsize What is the size of the population you would like to simulate? This is the population from which you will ultimately draw your samples. Default = 100,000. Accepts any positive whole number
#' @param iterations How many times you would like to estimate your models in random samples drawn from your population? One model will be run in each random sample. Default = 10,000. Accepts any whole number greater than 0.
#' @param alpha Set your alpha level. This is the threshold below which p-values will be considered significant. Default = 0.05. Accepts any number.
#' @param seed Set a seed to make your results reproducible. Default = 1. Accepts any number.
#' @param n Set the size of each sample to be drawn from the population. This is the sample size for which you are estimating statistical power. In other words, setting n to equal 100 will estimate statistical power at n = 100. Accepts any positive whole number smaller than your population.
#' @param d1...d10 The effect size, expressed as Cohen's d, for each model. You should specify a number of effect sizes equal to the number you specified in the "models" argument. That is, if you set models to equal 4, you should supply values for d1, d2, d3, and d4. You should always specify these in order, beginning with d1, and not skipping any. You must specify effect sizes for at least d1 and d2. Accepts any number.
#'
#' @value A dataframe containing a power estiamte, expressed as a decimal, for each of the models individually, and for all the models simultaneously.
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
#' simusim.bivar(n = 150, d1 = .2, d2 = .45)
#'
#' # Another example, customizing additional parameters.
#'
#' simusim.bivar(n = 300, d1 = .6, d2 = .43, d3 = .21, models = 3,
#'      popsize = 500000, alpha = .01, seed = 123)
#'
#' @export

simusim.bivar <- function(n, d1, d2, 
                          d3 = 0, d4 = 0, d5 = 0, d6 = 0, d7 = 0, d8 = 0, d9 = 0, d10 = 0,
                          models = 2, popsize = 100000, iterations = 10000, alpha = .05, seed = 1) {
  
  # Throw a warning if the user has specified the wrong number of models #
  dummy_es <- 0
  
  specified_models <- table(c(d1, d2, d3, d4, d5, d6, d7, d8, d9, d10, dummy_es))
  
  if(specified_models["0"] != (11 - models)) {
    stop("You have specified the wrong number of models. Make sure to set the models argument to equal the number of models for which you have specified effect sizes. Additionally, ensure to specify effect sizes for models in order, beginning with d1, and not skipping any along the way. ")
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
  
  # Throw a warning if the sample size is 0 or negative #
  if(popsize < 1 | round(popsize, 0) != popsize) {
    stop("You have specified an invalid population size. Please specify a whole number greater than 0.")
  }

  # Simulate the population #
  cortable <- matrix(c(1,d1/(sqrt(d1**2 + 4)),0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                       d1/(sqrt(d1**2 + 4)),1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                       0,0,1,d2/(sqrt(d2**2 + 4)),0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                       0,0,d2/(sqrt(d2**2 + 4)),1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                       0,0,0,0,1,d3/(sqrt(d3**2 + 4)),0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                       0,0,0,0,d3/(sqrt(d3**2 + 4)),1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                       0,0,0,0,0,0,1,d4/(sqrt(d4**2 + 4)),0,0,0,0,0,0,0,0,0,0,0,0,
                       0,0,0,0,0,0,d4/(sqrt(d4**2 + 4)),1,0,0,0,0,0,0,0,0,0,0,0,0,
                       0,0,0,0,0,0,0,0,1,d5/(sqrt(d5**2 + 4)),0,0,0,0,0,0,0,0,0,0,
                       0,0,0,0,0,0,0,0,d5/(sqrt(d5**2 + 4)),1,0,0,0,0,0,0,0,0,0,0,
                       0,0,0,0,0,0,0,0,0,0,1,d6/(sqrt(d6**2 + 4)),0,0,0,0,0,0,0,0,
                       0,0,0,0,0,0,0,0,0,0,d6/(sqrt(d6**2 + 4)),1,0,0,0,0,0,0,0,0,
                       0,0,0,0,0,0,0,0,0,0,0,0,1,d7/(sqrt(d7**2 + 4)),0,0,0,0,0,0,
                       0,0,0,0,0,0,0,0,0,0,0,0,d7/(sqrt(d7**2 + 4)),1,0,0,0,0,0,0,
                       0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,d8/(sqrt(d8**2 + 4)),0,0,0,0,
                       0,0,0,0,0,0,0,0,0,0,0,0,0,0,d8/(sqrt(d8**2 + 4)),1,0,0,0,0,
                       0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,d9/(sqrt(d9**2 + 4)),0,0,
                       0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,d9/(sqrt(d9**2 + 4)),1,0,0,
                       0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,d10/(sqrt(d10**2 + 4))), 20, 20)
  
  set.seed(seed)
  suppressWarnings(invisible(capture.output(dist <- SimCorrMix::corrvar(n = popsize,
                                                       k_cont = 20,
                                                       method = "Fleishman",
                                                       means = rep(0, times = 20),
                                                       vars = rep(1, times = 20),
                                                       skews = rep(0, times = 20),
                                                       skurts = rep(0, times = 20),
                                                       rho = cortable)$Y_cont)))
  
  dist <- data.frame(dist)
  
  names(dist) <- c("iv1", "dv1",
                   "iv2", "dv2",
                   "iv3", "dv3",
                   "iv4", "dv4",
                   "iv5", "dv5",
                   "iv6", "dv6",
                   "iv7", "dv7",
                   "iv8", "dv8",
                   "iv9", "dv9",
                   "iv10", "dv10")
  
  # Let the user know it's working #
  message(paste("Running", iterations, "simulations. This may take a minute.", sep = " "))
  
  # Run the models i times #
  result <- vector()
  
  set.seed(seed + 1)
  
  for(i in 1:iterations) {
    
    sample <- dist[sample(nrow(dist), size = n),]
    
    if (models == 2) {
      model1 <- data.frame(summary(lm(dv1 ~ iv1, data = sample))$coefficients)
      model2 <- data.frame(summary(lm(dv2 ~ iv2, data = sample))$coefficients)
      
      d1_result <- model1$Pr...t..[2] < alpha
      d2_result <- model2$Pr...t..[2] < alpha

      round_results <- data.frame(d1_result, d2_result)
      result <- rbind(result, 
                      round_results)
    } else if (models == 3) {
      model1 <- data.frame(summary(lm(dv1 ~ iv1, data = sample))$coefficients)
      model2 <- data.frame(summary(lm(dv2 ~ iv2, data = sample))$coefficients)
      model3 <- data.frame(summary(lm(dv3 ~ iv3, data = sample))$coefficients)
      
      d1_result <- model1$Pr...t..[2] < alpha
      d2_result <- model2$Pr...t..[2] < alpha
      d3_result <- model3$Pr...t..[2] < alpha
      
      round_results <- data.frame(d1_result, d2_result, d3_result)
      result <- rbind(result, 
                      round_results)
    } else if (models == 4) {
      model1 <- data.frame(summary(lm(dv1 ~ iv1, data = sample))$coefficients)
      model2 <- data.frame(summary(lm(dv2 ~ iv2, data = sample))$coefficients)
      model3 <- data.frame(summary(lm(dv3 ~ iv3, data = sample))$coefficients)
      model4 <- data.frame(summary(lm(dv4 ~ iv4, data = sample))$coefficients)
      
      d1_result <- model1$Pr...t..[2] < alpha
      d2_result <- model2$Pr...t..[2] < alpha
      d3_result <- model3$Pr...t..[2] < alpha
      d4_result <- model4$Pr...t..[2] < alpha
      
      round_results <- data.frame(d1_result, d2_result, d3_result, d4_result)
      result <- rbind(result, 
                      round_results)
    } else if (models == 5) {
      model1 <- data.frame(summary(lm(dv1 ~ iv1, data = sample))$coefficients)
      model2 <- data.frame(summary(lm(dv2 ~ iv2, data = sample))$coefficients)
      model3 <- data.frame(summary(lm(dv3 ~ iv3, data = sample))$coefficients)
      model4 <- data.frame(summary(lm(dv4 ~ iv4, data = sample))$coefficients)
      model5 <- data.frame(summary(lm(dv5 ~ iv5, data = sample))$coefficients)
      
      d1_result <- model1$Pr...t..[2] < alpha
      d2_result <- model2$Pr...t..[2] < alpha
      d3_result <- model3$Pr...t..[2] < alpha
      d4_result <- model4$Pr...t..[2] < alpha
      d5_result <- model5$Pr...t..[2] < alpha
      
      round_results <- data.frame(d1_result, d2_result, d3_result, d4_result, d5_result)
      result <- rbind(result, 
                      round_results)
    } else if (models == 6) {
      model1 <- data.frame(summary(lm(dv1 ~ iv1, data = sample))$coefficients)
      model2 <- data.frame(summary(lm(dv2 ~ iv2, data = sample))$coefficients)
      model3 <- data.frame(summary(lm(dv3 ~ iv3, data = sample))$coefficients)
      model4 <- data.frame(summary(lm(dv4 ~ iv4, data = sample))$coefficients)
      model5 <- data.frame(summary(lm(dv5 ~ iv5, data = sample))$coefficients)
      model6 <- data.frame(summary(lm(dv6 ~ iv6, data = sample))$coefficients)
      
      d1_result <- model1$Pr...t..[2] < alpha
      d2_result <- model2$Pr...t..[2] < alpha
      d3_result <- model3$Pr...t..[2] < alpha
      d4_result <- model4$Pr...t..[2] < alpha
      d5_result <- model5$Pr...t..[2] < alpha
      d6_result <- model6$Pr...t..[2] < alpha
      
      round_results <- data.frame(d1_result, d2_result, d3_result, d4_result, d5_result, d6_result)
      result <- rbind(result, 
                      round_results)
    } else if (models == 7) {
      model1 <- data.frame(summary(lm(dv1 ~ iv1, data = sample))$coefficients)
      model2 <- data.frame(summary(lm(dv2 ~ iv2, data = sample))$coefficients)
      model3 <- data.frame(summary(lm(dv3 ~ iv3, data = sample))$coefficients)
      model4 <- data.frame(summary(lm(dv4 ~ iv4, data = sample))$coefficients)
      model5 <- data.frame(summary(lm(dv5 ~ iv5, data = sample))$coefficients)
      model6 <- data.frame(summary(lm(dv6 ~ iv6, data = sample))$coefficients)
      model7 <- data.frame(summary(lm(dv7 ~ iv7, data = sample))$coefficients)
      
      d1_result <- model1$Pr...t..[2] < alpha
      d2_result <- model2$Pr...t..[2] < alpha
      d3_result <- model3$Pr...t..[2] < alpha
      d4_result <- model4$Pr...t..[2] < alpha
      d5_result <- model5$Pr...t..[2] < alpha
      d6_result <- model6$Pr...t..[2] < alpha
      d7_result <- model7$Pr...t..[2] < alpha
      
      round_results <- data.frame(d1_result, d2_result, d3_result, d4_result, d5_result, d6_result, d7_result)
      result <- rbind(result, 
                      round_results)
    } else if (models == 8) {
      model1 <- data.frame(summary(lm(dv1 ~ iv1, data = sample))$coefficients)
      model2 <- data.frame(summary(lm(dv2 ~ iv2, data = sample))$coefficients)
      model3 <- data.frame(summary(lm(dv3 ~ iv3, data = sample))$coefficients)
      model4 <- data.frame(summary(lm(dv4 ~ iv4, data = sample))$coefficients)
      model5 <- data.frame(summary(lm(dv5 ~ iv5, data = sample))$coefficients)
      model6 <- data.frame(summary(lm(dv6 ~ iv6, data = sample))$coefficients)
      model7 <- data.frame(summary(lm(dv7 ~ iv7, data = sample))$coefficients)
      model8 <- data.frame(summary(lm(dv8 ~ iv8, data = sample))$coefficients)
      
      d1_result <- model1$Pr...t..[2] < alpha
      d2_result <- model2$Pr...t..[2] < alpha
      d3_result <- model3$Pr...t..[2] < alpha
      d4_result <- model4$Pr...t..[2] < alpha
      d5_result <- model5$Pr...t..[2] < alpha
      d6_result <- model6$Pr...t..[2] < alpha
      d7_result <- model7$Pr...t..[2] < alpha
      d8_result <- model8$Pr...t..[2] < alpha
      
      round_results <- data.frame(d1_result, d2_result, d3_result, d4_result, d5_result, d6_result, d7_result, d8_result)
      result <- rbind(result, 
                      round_results)
    } else if (models == 9) {
      model1 <- data.frame(summary(lm(dv1 ~ iv1, data = sample))$coefficients)
      model2 <- data.frame(summary(lm(dv2 ~ iv2, data = sample))$coefficients)
      model3 <- data.frame(summary(lm(dv3 ~ iv3, data = sample))$coefficients)
      model4 <- data.frame(summary(lm(dv4 ~ iv4, data = sample))$coefficients)
      model5 <- data.frame(summary(lm(dv5 ~ iv5, data = sample))$coefficients)
      model6 <- data.frame(summary(lm(dv6 ~ iv6, data = sample))$coefficients)
      model7 <- data.frame(summary(lm(dv7 ~ iv7, data = sample))$coefficients)
      model8 <- data.frame(summary(lm(dv8 ~ iv8, data = sample))$coefficients)
      model9 <- data.frame(summary(lm(dv9 ~ iv9, data = sample))$coefficients)
      
      d1_result <- model1$Pr...t..[2] < alpha
      d2_result <- model2$Pr...t..[2] < alpha
      d3_result <- model3$Pr...t..[2] < alpha
      d4_result <- model4$Pr...t..[2] < alpha
      d5_result <- model5$Pr...t..[2] < alpha
      d6_result <- model6$Pr...t..[2] < alpha
      d7_result <- model7$Pr...t..[2] < alpha
      d8_result <- model8$Pr...t..[2] < alpha
      d9_result <- model9$Pr...t..[2] < alpha
      
      round_results <- data.frame(d1_result, d2_result, d3_result, d4_result, d5_result, d6_result, d7_result, d8_result, d9_result)
      result <- rbind(result, 
                      round_results)
    } else if (models == 10) {
      model1 <- data.frame(summary(lm(dv1 ~ iv1, data = sample))$coefficients)
      model2 <- data.frame(summary(lm(dv2 ~ iv2, data = sample))$coefficients)
      model3 <- data.frame(summary(lm(dv3 ~ iv3, data = sample))$coefficients)
      model4 <- data.frame(summary(lm(dv4 ~ iv4, data = sample))$coefficients)
      model5 <- data.frame(summary(lm(dv5 ~ iv5, data = sample))$coefficients)
      model6 <- data.frame(summary(lm(dv6 ~ iv6, data = sample))$coefficients)
      model7 <- data.frame(summary(lm(dv7 ~ iv7, data = sample))$coefficients)
      model8 <- data.frame(summary(lm(dv8 ~ iv8, data = sample))$coefficients)
      model9 <- data.frame(summary(lm(dv9 ~ iv9, data = sample))$coefficients)
      model10 <- data.frame(summary(lm(dv10 ~ iv10, data = sample))$coefficients)
      
      d1_result <- model1$Pr...t..[2] < alpha
      d2_result <- model2$Pr...t..[2] < alpha
      d3_result <- model3$Pr...t..[2] < alpha
      d4_result <- model4$Pr...t..[2] < alpha
      d5_result <- model5$Pr...t..[2] < alpha
      d6_result <- model6$Pr...t..[2] < alpha
      d7_result <- model7$Pr...t..[2] < alpha
      d8_result <- model8$Pr...t..[2] < alpha
      d9_result <- model9$Pr...t..[2] < alpha
      d10_result <- model10$Pr...t..[2] < alpha
      
      round_results <- data.frame(d1_result, d2_result, d3_result, d4_result, d5_result, d6_result, d7_result, d8_result, d9_result, d10_result)
      result <- rbind(result, 
                      round_results)
    }
    
  }
  
  # Create variable for calculating simultaneous power #
  result$simultaneous <- FALSE
  
  if (models == 2) {
    result$simultaneous[result$d1_result == TRUE & result$d2_result == TRUE] <- TRUE
  } else if (models == 3) {
    result$simultaneous[result$d1_result == TRUE & result$d2_result == TRUE & result$d3_result == TRUE] <- TRUE
  } else if (models == 4) {
    result$simultaneous[result$d1_result == TRUE & result$d2_result == TRUE & result$d3_result == TRUE & result$d4_result == TRUE] <- TRUE
  } else if (models == 5) {
    result$simultaneous[result$d1_result == TRUE & result$d2_result == TRUE & result$d3_result == TRUE & result$d4_result == TRUE & result$d5_result == TRUE] <- TRUE
  } else if (models == 6) {
    result$simultaneous[result$d1_result == TRUE & result$d2_result == TRUE & result$d3_result == TRUE & result$d4_result == TRUE & result$d5_result == TRUE &
                          result$d6_result == TRUE] <- TRUE
  } else if (models == 7) {
    result$simultaneous[result$d1_result == TRUE & result$d2_result == TRUE & result$d3_result == TRUE & result$d4_result == TRUE & result$d5_result == TRUE &
                          result$d6_result == TRUE & result$d7_result == TRUE] <- TRUE
  } else if (models == 8) {
    result$simultaneous[result$d1_result == TRUE & result$d2_result == TRUE & result$d3_result == TRUE & result$d4_result == TRUE & result$d5_result == TRUE &
                          result$d6_result == TRUE & result$d7_result == TRUE & result$d8_result == TRUE] <- TRUE
  } else if (models == 9) {
    result$simultaneous[result$d1_result == TRUE & result$d2_result == TRUE & result$d3_result == TRUE & result$d4_result == TRUE & result$d5_result == TRUE &
                          result$d6_result == TRUE & result$d7_result == TRUE & result$d8_result == TRUE & result$d9_result == TRUE] <- TRUE
  } else if (models == 10) {
    result$simultaneous[result$d1_result == TRUE & result$d2_result == TRUE & result$d3_result == TRUE & result$d4_result == TRUE & result$d5_result == TRUE &
                          result$d6_result == TRUE & result$d7_result == TRUE & result$d8_result == TRUE & result$d9_result == TRUE & result$d10_result == TRUE] <- TRUE
  }
  
  result[result == "TRUE"] <- 1
  
  # Calculate power #
  simultaneous_power <- mean(result$simultaneous)
  model1_power <- mean(result$d1_result)
  model2_power <- mean(result$d2_result)
  if (models > 2) {
    model3_power <- mean(result$d3_result)
  }
  if (models > 3) {
    model4_power <- mean(result$d4_result)
  }
  if (models > 4) {
    model5_power <- mean(result$d5_result)
  }
  if (models > 5) {
    model6_power <- mean(result$d6_result)
  }
  if (models > 6) {
    model7_power <- mean(result$d7_result)
  }
  if (models > 7) {
    model8_power <- mean(result$d8_result)
  }
  if (models > 8) {
    model9_power <- mean(result$d9_result)
  }
  if (models > 9) {
    model10_power <- mean(result$d10_result)
  }
  
  # Summarize results
  
  if (models == 2) {
    power <- data.frame(model1_power, model2_power, simultaneous_power)
  } else if ( models == 3) {
    power <- data.frame(model1_power, model2_power, model3_power, simultaneous_power)
  } else if ( models == 4) {
    power <- data.frame(model1_power, model2_power, model3_power, model4_power, simultaneous_power)
  } else if ( models == 5) {
    power <- data.frame(model1_power, model2_power, model3_power, model4_power, model5_power, simultaneous_power)
  } else if ( models == 6) {
    power <- data.frame(model1_power, model2_power, model3_power, model4_power, model5_power, model6_power, simultaneous_power)
  } else if ( models == 7) {
    power <- data.frame(model1_power, model2_power, model3_power, model4_power, model5_power, model6_power, model7_power, simultaneous_power)
  } else if ( models == 8) {
    power <- data.frame(model1_power, model2_power, model3_power, model4_power, model5_power, model6_power, model7_power, model8_power, simultaneous_power)
  } else if ( models == 9) {
    power <- data.frame(model1_power, model2_power, model3_power, model4_power, model5_power, model6_power, model7_power, model8_power, model9_power, simultaneous_power)
  } else if ( models == 10) {
    power <- data.frame(model1_power, model2_power, model3_power, model4_power, model5_power, model6_power, model7_power, model8_power, model9_power, model10_power, simultaneous_power)
  }
  
  return(power)
}
