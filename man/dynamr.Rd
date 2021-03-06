% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dynamr.R
\name{dynamr}
\alias{dynamr}
\title{A Permutation-Based Changepoint Technique for Monitoring Effect Sizes}
\usage{
dynamr(dat, time_var, formula, N = 2000, window_size,
  family = c("gaussian", "binomial", "poisson"))
}
\arguments{
\item{dat}{Data frame with panel data of interest}

\item{time_var}{Name of variable denoting time period for any observation,
must be a character}

\item{formula}{Formula for glm object}

\item{N}{Number of coefficient samples extracted for calculating stable
statistical behavior, defaults to 2000}

\item{window_size}{Number of time periods per model}

\item{family}{A description of the error distribution to be used for the
fitted glm. Currently, there are three options: "gaussian" - multivariate
 linear regression; "binomial" - multivariate logistic regression; "poisson"
 -- multivariate Poisson regression}
}
\description{
\code{dynamr} is used to detect where effect values change temporally when
estimating generalized linear models with panel data.
}
\examples{
library(ISLR)
library(dynamr)

data("Weekly")

stock_dynam <- dynamr(
 dat = Weekly,
 time_var = "Year",
 formula = Today ~ Lag1 + Lag2,
 window_size = 1,
 family = "gaussian",
 N = 5000
)

stock_dynam

}
