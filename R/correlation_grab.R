##' This function calculates the resulting serial correlation between grab samples each having \code{r} primary increments with original serial correlation \code{d}.
##' @title Serial correlation between grab samples
##' @param r number of primary increments in a grab sample or grab sample size
##' @param p limiting fraction or proportion of contaminated increments
##' @param d serial correlation of contamination between the primary increments
##' @return Serial correlation between grab samples
##' @details The serial correlation between blocks (grab samples) is given by \eqn{d_g} as
##' \deqn{d_g = [dp(1-p(1-d))^{r-1}]/p_d}
##' where \eqn{p_d} is the probability of detection in any of the block (grab sample) which is calculated by using \link{prob_detect_single_grab}.
##' @seealso \link{prob_detect_single_grab}
##' @examples
##' r <-  25
##' p <-  0.005
##' d <-  0.99
##' correlation_grab(r, p, d)
##' @export
## we have used the notation as d_g in the paper Quantitative risk assessment for grab sampling inspection of powdered products
correlation_grab <- function(r, p, d) {
  result <- (d * p * (1 - p * (1 - d))^(r - 1))/prob_detect_single_grab(r, p, d)
  return(result)
  }


