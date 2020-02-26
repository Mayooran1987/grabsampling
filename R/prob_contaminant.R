##' This function calculates the probability of exactly  \code{l}  contaminated samples out of \code{t} selected grab samples for given gram sample size \code{r} and serial correlation \code{d} at the process contamination level \code{p} for a production length of \code{N}.
##' @title Probability of contaminated sample
##' @param l number of contaminated in \code{t} selected samples
##' @param r number of primary increments in a grab sample or grab sample size
##' @param t number of grab samples
##' @param d serial correlation of contamination between the primary increments
##' @param p limiting fraction or proportion of contaminated increments
##' @param N length of the production
##' @param method what sampling method we have applied such as \code{'systematic'} or \code{'random'} selection methods
##' @return Probability of contaminated
##' @details Let \eqn{S_t} be the number of contaminated samples and \eqn{S_t=\sum X_t} where \eqn{X_t=1} or \eqn{0} depending on the presence or absence of contamination, then \eqn{P(S_t=l)} formula given in \href{https://doi.org/10.2307/1427041}{Bhat and Lal (1988)}, also we can use following recurrence relation formula,
##' \deqn{P(S_t=l)=P(X_t=1;S_{t-1}=l-1) + P(X_t=0;S_{t-1}=l)} which is given in \href{https://onlinelibrary.wiley.com/doi/abs/10.1002/nav.1028}{Vellaisamy and Sankar (2001)}. Both methods will be produced the same results.
##' For this package development, we directly applied formula which is from \href{https://doi.org/10.2307/1427041}{Bhat and Lal (1988)}.
##' @seealso \link{prob_detect_single_grab}, \link{correlation_grab}
##' @references
##' \itemize{
##' \item  Bhat, U., & Lal, R. (1988). Number of successes in Markov trials. Advances in Applied Probability, 20(3), \href{https://doi.org/10.2307/1427041}{677-680}.
##' \item  Vellaisamy, P., Sankar, S., (2001). Sequential and systematic sampling plans for the Markov-dependent production process. Naval Research Logistics 48, \href{https://onlinelibrary.wiley.com/doi/abs/10.1002/nav.1028}{451-467}.
##' }
##' @examples
##'   l <-  1
##'   r <-  25
##'   t <-  30
##'   d <-  0.99
##'   p <-  0.005
##'   N <-  1e9
##'   method <- 'systematic'
##'   prob_contaminant(l, r, t, d, p, N, method)
##' @export
prob_contaminant <- function(l, r, t, d, p, N, method) {
    p_d <- prob_detect_single_grab(r, p, d)
    if (method == "random") {
        s_t <- stats::dbinom(l, t, p_d)
        } else if (method == "systematic") {
            d_g <- correlation_grab(r, p, d)
            k <- ceiling(N/(r * t))
            sum1 <- 0
            sum2 <- 0
            if (l == 0) {
                sum1 <- (1 - p_d) * (1 - (p_d * (1 - d_g^k)))^(t - 1)
                sum2 <- 0
                sum3 <- 0
                } else if (l == 1) {
                    sum1 <- (1 - p_d) * ((1 - (p_d * (1 - d_g^k)))^(t - 3)) * (p_d * (1 - d_g^k)) * (((1 - p_d) * (1 - d_g^k) * choose(t - 2, 1)) + (1 - (p_d *(1 - d_g^k))) * choose(t - 2, 0))
                    sum2 <- 0
                    sum3 <- p_d * ((1 - p_d) * (1 - d_g^k)) * ((1 - (p_d * (1 - d_g^k)))^(t - 2))
                    } else {
                        for (j in min(1, l):l) {
                            sum1 <- sum1 + (1 - p_d) * ((choose(l - 1, j - 1)) * ((1 - (p_d * (1 - d_g^k)))^(t - l - j - 1)) * ((1 - ((1 - p_d) * (1 - d_g^k)))^(l -j)) * (p_d * (1 - d_g^k))^j * (((1 - p_d) * (1 - d_g^k))^(j - 1)) * (((1 - p_d) * (1 - d_g^k)) * (ifelse(t - l  < t , 1, 0)) *
                                                                                                                   choose(t - l - 1, j) + (1 - (p_d * (1 - d_g^k))) *(ifelse(l == 0, 1, 0)) * (ifelse(l + j - 1 < t, 1, 0)) * choose(t - l - 1, j - 1)))
                            }
                        for (i in 1:(l - 1)) {
                                sum2 <- sum2 + p_d * (choose(l - 1, i)) * ((1 - (p_d * (1 - d_g^k)))^(t - l - i - 1)) * ((1 - ((1 - p_d) * (1 - d_g^k)))^(l - 1 - i)) *(p_d * (1 - d_g^k))^i * ((1 - p_d) * (1 - d_g^k))^i * (((1 - p_d) * (1 - d_g^k)) * (ifelse(l + i < t , 1, 0)) * choose(t - l -1, i) + (1 - (p_d * (1 - d_g^k))) *(ifelse(l + i - 2 < t - 1, 1, 0)) * choose(t - l - 1, i - 1))
                            }
                        }
            sum3 <- p_d * (ifelse(l == 0, 1, 0)) * ((1 - ((1 - p_d) * (1 - d_g^k)))^(l - 1)) * ((ifelse(l == t , 1, 0)) +  (ifelse(l == 0, 1, 0) * ifelse(t - 1 == 0, 1, 0) * (1 - ifelse(t - 1 < l, 1, 0)))  * ((1 - p_d) *(1 - d_g^k)) * ((1 - (p_d * (1 - d_g^k)))^(t - 1 - l)))
            s_t <- sum1 + sum2 + sum3
        }
    else {
        warning ("please choose one of the given sampling method with case sensitive such as 'random' or 'systematic'")
    }
    return(s_t)
    }
## We have derived transition matrix for sequential sampling method, so we can include sequential sampling method too if we want.

