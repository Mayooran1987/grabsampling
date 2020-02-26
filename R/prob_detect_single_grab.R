##' This function calculates the probability of detection in a single grab sample comprising of \code{r} primary increments for given serial correlation \code{d}.
##' @title Probability of detection in a single grab sample
##' @param r number of primary increments in a grab sample or grab sample size
##' @param p limiting fraction or proportion of contaminated increments
##' @param d serial correlation of contamination between the primary increments
##' @details The probability of detection in any of the grab sample is given by \eqn{p_d} as
##' \deqn{p_d = 1-(1-p)(1-p(1-d))^{r-1}}
##' @return Probability of detection in a grab sample
##' @examples
##'    r <-  25
##'    p <-  0.005
##'    d <-  0.99
##'    prob_detect_single_grab(r, p, d)
##' @export
## we have used the notation as p_d in the paper Quantitative risk assessment for grab sampling inspection of powdered products
prob_detect_single_grab <- function(r, p, d) {
  result <- 1 - (1 - p) * (1 - p * (1 - d))^(r - 1)
  return(result)
  }

