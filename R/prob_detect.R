##' This function gives the detection probability for  \code{t} grab samples and given acceptance number under systematic or random sampling methods. This function is also used to calculate the detection probability for primary increments selection by setting the number of primary increments as one.
##' @title Probability of detection under the grab sampling method
##' @param c acceptance number
##' @param r number of primary increments in a grab sample or grab sample size
##' @param t number of grab samples
##' @param d serial correlation of contamination between the primary increments
##' @param p limiting fraction or proportion of contaminated increments
##' @param N length of the production
##' @param method what sampling method we have applied such as \code{'systematic'} or \code{'random'} selection methods
##' @return Probability of detection in all seleceted grab samples
##' @details The detection probability of entire selected grab samples is given by,
##' \deqn{P_D=1-[P(S_t=0)+P(S_t=1)+\cdots +P(S_t=c)]}
##' @seealso \link{prob_contaminant}
# @references None
##' @examples
##'   c <-  1
##'   r <-  25
##'   t <-  30
##'   d <-  0.99
##'   p <-  0.005
##'   N <-  1e9
##'   method <- 'systematic'
##'   prob_detect(c, r, t, d, p, N, method)
##' @export
## we have used the notation as P_D in the paper Quantitative risk assessment for grab sampling inspection of powdered products
prob_detect <- function(c, r, t, d, p, N, method) {
    prob <- 0
    for (i in 0:c) {
        prob <- prob + prob_contaminant(i, r, t, d, p, N, method)
        }
    prob_detection <- 1 - prob
    return(prob_detection)
    }
