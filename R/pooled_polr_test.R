#' Test an orderd logistic or probit regression for the presence of time-varying effects
#'
#' \code{pooled_polr_test} is used to test the temporal effect homogeneity
#' of a ordered logistic or probit regression with panel data.
#'
#' @param formula An object of class "formula", which provides a symbolic
#' description of the model to be fitted. Ensure that names included here are
#' contained with the data frame entries of the Panel.data object.
#' @param Panel.data A list of data frames containing the predictors and
#' response variables listed in the model formula. List entries are assumed
#' to be ordered according to time.
#' @param N The number of bootstrap samples to run for this test. Default is 100.
#' @param method Logistic, probit, log-log, or cauchit
#'
#' @return
#' @export
#'

pooled_polr_test <- function(
  formula,
  Panel.data,
  N = 100,
  method = c(
    "logistic",
    "probit",
    "loglog",
    "cloglog",
    "cauchit"
  )
){
  # Step 1 - calculate Gamma_0 = generalized likelihood ratio statistic
  GLR <- function(
    formula,
    Panel.data,
    method
  ){
    # Run pooled data fit first
    pooled.data <- Reduce("rbind", Panel.data)
    pooled.model <- polr(
      formula = formula,
      method = method,
      data = pooled.data
    )
    likelihood.pooled <- as.numeric(logLik(pooled.model))

    # Run varying coefficient model now
    T <- length(Panel.data)
    likelihood.vc <- 0
    for(t in 1:T){
      model.temp <- polr(
        formula = formula,
        method = method,
        data = Panel.data[[t]]
      )
      likelihood.vc <- likelihood.vc + as.numeric(logLik(model.temp))
    }
    return(Gamma <- likelihood.vc - likelihood.pooled)
  }
  # Calculate generalized likelihood ratio statistic for observed data
  Gamma0 <- GLR(formula = formula, Panel.data = Panel.data, method = method)

  # Calculate number of observations in each data set
  n.obs <- sapply(Panel.data, function(x){dim(x)[1]})

  # Calculate the pooled model so that we can simulate from it
  pooled.data <- Reduce("rbind", Panel.data)
  pooled.model <- polr(
    formula = formula,
    method = method,
    data = pooled.data
  )

  T <- length(Panel.data)
  Gamma <- rep(0, N)

  # Run Bootstrap procedure to calculate N GLR statistics
  for(n in 1:N){
    # create new y values according to the length of each data set
    Panel.new <- Panel.data
    response.name <- as.character(formula)[2]
    response.indx <- which(names(Panel.new[[1]]) == response.name)

    # simulate new values (according to family)
    y.new <- simulate(pooled.model, type = "response")
    # place values into each time point
    for(t in 1:T){
      if(t == 1){
        Panel.new[[t]][, response.indx] <- y.new[1:n.obs[t], ]
      }
      if(t > 1){
        Panel.new[[t]][, response.indx] <- y.new[
          (n.obs[t-1]+1):(n.obs[t-1] + n.obs[t]),
          ]
      }
    }
    Gamma[n] <- GLR(
      formula = formula,
      Panel.data = Panel.new,
      method = method
    )
    print(n)
  }
  p.value <- length(which(Gamma > Gamma0))/N
  return(list(
    Simulated.Gamma = Gamma,
    Observed.Gamma = Gamma0,
    p.value = p.value
  )
  )
}
