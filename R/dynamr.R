#' A Permutation-Based Changepoint Technique for Monitoring Effect Sizes
#'
#' \code{dynamr} is used to detect where effect values change temporally when
#' estimating generalized linear models with panel data.
#'
#' @param dat Data frame with panel data of interest
#' @param time_var Name of variable denoting time period for any observation,
#' must be a character
#' @param formula Formula for glm object
#' @param N Number of coefficient samples extracted for calculating stable
#' statistical behavior, defaults to 2000
#' @param window_size Number of time periods per model
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
#' library(dynamr)
#'
#' data("Weekly")
#'
#' stock_dynam <- dynamr(
#'  dat = Weekly,
#'  time_var = "Year",
#'  formula = Today ~ Lag1 + Lag2,
#'  window_size = 1,
#'  family = "gaussian",
#'  N = 5000
#' )
#'
#' stock_dynam
#'
dynamr <- function(
  dat,
  time_var,
  formula,
  N = 2000,
  window_size,
  family = c(
    "gaussian",
    "binomial",
    "poisson"
  )
) {
  require(dplyr)
  require(tibble)

  # Denote start and end times
  start_time <- min(dat[[time_var]], na.rm = TRUE)
  end_time <- max(dat[[time_var]], na.rm = TRUE)
  
  # Number of models to fit over time
  # + 2 works out to include all time periods
  num_models <- end_time - (start_time + window_size) + 2
  
  # Storage for coefficient estimates
  coefs <- matrix(
    NA,
    nrow = N,
    ncol = length(all.vars(formula[-2])) + 1 # -2 removes DV, +1 intercept
  )
  
  # Create dist of time-invariant coefficients
  # For each variable
  for (i in 1:N) {
    t <- sample(c(start_time:end_time), window_size)
    keep_obs <- which(dat[[time_var]] %in% t)
    mod <- glm( # model with randomly selected time periods
      formula,
      data = dat[keep_obs, ],
      family = family
    )
    ## Which to extract?
    coefs[i, ] <- coef(mod)
  }
  
  # Create ncol tolerance regions (shewhart)
  region <- 3 * apply(coefs, MARGIN = 2, sd)
  
  # 3 columns per var
  # extract variable name for column names
  changes <- tibble(
    "time" = as.integer(),
    "int_est" = as.numeric(),
    "int_sd" = as.numeric(),
    "pr_int_change" = as.numeric()
  )
  
  # Automate dynamic dynamic by formula terms
  for (i in 1:length(all.vars(formula[-2]))) {
    changes <- changes %>%
      add_column(!!
                   paste0(all.vars(formula[-2])[i], "_est", sep = ""),
      ) %>%
      add_column(!!
                   paste0(all.vars(formula[-2])[i], "_se", sep = ""),
      ) %>%
      add_column(!!
                   paste0("pr_", all.vars(formula[-2])[i], "_change", sep = ""),
      )
  }
  
  # Make columns are numeric vars
  changes <- changes %>%
    mutate(across(where(is.character), as.numeric))
  
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
    est <- summary(mod)$coefficients[,1:2]
    
    if (i == 1) {
      mu <- est[, 1] # coefficients
      mu_upper <- mu + region # upper bounds
      mu_lower <- mu - region # lower bounds
      
      # Append time -- start the window with middle time
      changes[i, 1] <- floor(median(seq(time_min, time_max)))
      
      # Fill in estimates df
      # Will work for any number of covariates
      for (j in 1:nrow(est)) {
        changes[i, 3 * j - 1] <- est[j, 1] # coefficient
        changes[i, 3 * j] <- est[j, 2] # sd
        changes[i, 3 * j + 1] <- 0 # first point can"t be change
      }
    }else{
      # append time
      changes[i, 1] <- floor(median(seq(time_min, time_max)))
      
      # estimates
      for (j in 1:nrow(est)) {
        changes[i, 3 * j - 1] <- est[j, 1] # coefficient
        changes[i, 3 * j] <- est[j, 2] # sd
        # pr(observed coef is larger than window"s upper bound)
        pr_above_upper <- 1 - pnorm(
          as.numeric(mu_upper[j]),
          mean = as.numeric(changes[i, 3 * j - 1]),
          sd = as.numeric(changes[i, 3 * j])
        )
        pr_below_lower <- pnorm(
          as.numeric(mu_lower[j]),
          mean = as.numeric(changes[i, 3 * j - 1]),
          sd = as.numeric(changes[i, 3 * j])
        )
        # Record prob of change
        ifelse(
          pr_above_upper > pr_below_lower,
          changes[i, 3 * j + 1] <- pr_above_upper,
          changes[i, 3 * j + 1] <- pr_below_lower
        )
        # update bounds if so
        if (pr_above_upper >= 0.5 | pr_below_lower >= 0.5) {
          mu[j] <- as.numeric(changes[i, 3 * j - 1]) # is a change
          mu_upper[j] <- mu[j] + region[j]
          mu_lower[j] <- mu[j] - region[j]
        }
      }
    }
  }
  return(changes)
}
