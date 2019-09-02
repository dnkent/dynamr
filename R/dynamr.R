#' A Randomization Approach to Dynamic Analysis of Panel Data
#'
#' \code{dynamr} is used to detect where effect values change temporally when
#' estimating generalized linear models with panel data.
#'
#' @param dat Data frame with panel data of interest
#' @param time_var Name of variable denoting time period for any observation,
#' must be a character
#' @param formula Formula for glm object
#' @param num_models Number of time-specific models to be ran
#' @param N Number of coefficient samples extracted for calculating stable
#' statistical behavior, defaults to 2000
#' @param window_size Number of time periods per model
#' @param start_time First time period value
#' @param end_time Last time period value
#' @param change_var Variable being checked for time-varying behavior,
#' must be a character
#' @param family A description of the error distribution to be used for the
#' fitted glm. Currently, there are three options: "gaussian" - multivariate
#'  linear regression; "binomial" - multivariate logistic regression; "poisson"
#'  -- multivariate Poisson regression
#'
#' @return
#' @export
#'
#' @examples
#' library(ISLR)
#' library(dplyr)
#'
#' data("Weekly")
#'
#' stock_dynam <- dynamr(
#'  dat = Weekly,
#'  time_var = "Year",
#'  formula = Today ~ Lag1 + Lag2,
#'  num_models = 21,
#'  window_size = 1,
#'  start_time = 1990,
#'  end_time = 2010,
#'  change_var = "Lag2",
#'  family = "gaussian",
#'  N = 5000
#' )
#'
#' stock_dynam
#'
#' stock_dynam %>% filter(change == 1)

dynamr <- function(
  dat,
  time_var,
  formula,
  num_models,
  N = 2000,
  window_size,
  start_time,
  end_time,
  change_var,
  family = c(
    "gaussian",
    "binomial",
    "poisson"
  )
){
  # Coefficient vector
  coefs <- rep(NA, N)

  # Create dist of time-invariant coefficients
  for(i in 1:N){
    t <- sample(c(start_time:end_time), window_size)
    keep_obs <- which(dat[[time_var]] %in% t)
    mod <- glm( # model with randomly selected time periods
      formula,
      data = dat[keep_obs, ],
      family = family
    )
    ## Which to extract?
    coef_num <- which(names(coef(mod)) == change_var)
    coefs[i] <- coef(mod)[coef_num]
  }

  # Tolerance region
  region <- 3 * sd(coefs)

  # Change locations
  changes <- tibble(
    "Time" = as.integer(),
    "Estimate" = as.numeric(),
    "std.err" = as.numeric(),
    "change" = as.integer()
  )

  # Extract coefficients
  for(i in 1:num_models){
    time_min <- (start_time + i) - 1
    time_max <- start_time + window_size + i - 1
    temp_obs <- which(
      dat[[time_var]] >= time_min &
        dat[[time_var]] < time_max
    )
    temp <- dat[temp_obs, ]
    mod <- glm(
      formula,
      data = temp,
      family = family
    )
    ## Which to extract?
    coef_num <- which(names(coef(mod)) == change_var)
    est <- summary(mod)$coefficients[coef_num, 1:2]

    if(i == 1){ ## Establish starting region
      mu <- est[1]
      mu_upper <- mu + region
      mu_lower <- mu - region
      # append time -- start  the window
      changes[i, 1] <- median(seq(time_min, time_max))
      ## because using a window, use middle time
      changes[i, 2] <- mu # append estimate
      changes[i, 3] <- est[2] # append confint
      changes[i, 4] <- 0 # append changes
    }else{
      coef_mw <- est[1] ## moving window coefficient
      # append time
      changes[i, 1] <- median(seq(time_min, time_max))
      ## because using a window, use middle time
      changes[i, 2] <- coef_mw # append coefficient
      changes[i, 3] <- est[2] # append confint
      if(coef_mw > mu_upper | coef_mw < mu_lower){
        changes[i, 4] <- 1 # append change
        mu <- coef_mw
        mu_upper <- mu + region # set new region
        mu_lower <- mu - region # set new region
      }else{
        changes[i, 4] <- 0 # append change
      }
    }
  }
  return(changes)
}
