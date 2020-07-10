## dynamr

An R package for the dynamic analysis of panel data, `dynamr` serves two purposes. First, the `pooled_glm_test()` function can be used by a researcher to test whether or not the assumption of stable effects holds when estimating a generalized linear model with panel data. Second, the `dynamr()` function detects where effect changes occur. 

The reference for this method is **A Permutation-Based Changepoint Technique for Monitoring Effect Sizes**, co-authored with [James D. Wilson](http://jdwilson-statistics.com/) and [Skyler Cranmer](http://www.skylercranmer.net/). The manuscript is currently under review but available upon request. 

## Installation

To install `dynamr` use the following command and make sure to have `devtools` installed.

```
devtools::install_github("dnkent/dynamr")
```

## Description

This package contains two primary functions which are briefly described below. For any function named ```function```, type ```?function``` in R to get full documentation.

- `dynamr()`: detect at which time points coefficients in a generalized linear model change. This function produces a tibble with the following output: time points, coefficient estimates for each time point, standard errors for each coefficient estimate, and the probability that a changepoint occurs for each coefficient estimate at that time period.

```
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
```

- `pooled_glm_test()`: test whether the assumption temporal effect homogeneity holds for coefficients of interest in a generalized linear model using panel data. 

```
library(ISLR)
library(dynamr)

data("Weekly")

stock_list <- list()

for(i in 1:21){
  time <- 1989 + i
  stock_list[[i]] <- dplyr::filter(Weekly, Year == time)
}

stock.p.value <- pooled_glm_test(
  formula = Today ~ Lag1 + Lag2,
  Panel.data = stock_list,
  N = 300,
  family = "gaussian"
)

stock.p.value$p.value
```

## Questions or Bugs

Please send any comments, bugs, or questions to the developer Daniel Kent at kent.249@osu.edu. 
