##' This function calculates the overall probability of acceptance for given microbiological distribution such as lognormal.
##' @title Probability of acceptance for grab sampling scheme
##' @param c acceptance number
##' @param r number of primary increments in a grab sample or grab sample size
##' @param t number of grab samples
##' @param mu location parameter (mean log) of the Lognormal and Poisson-lognormal distributions on the log10 scale
##' @param distribution what suitable microbiological distribution we have used such as  \code{'Poisson gamma'} or \code{'Lognormal'}or \code{'Poisson lognormal'}
##' @param K dispersion parameter of the Poisson gamma distribution (default value 0.25)
##' @param m microbiological limit with default value zero, generally expressed as number of microorganisms in specific sample weight
##' @param sd standard deviation of the lognormal and Poisson-lognormal distributions on the log10 scale (default value 0.8)
##' @return Probability of acceptance
##' @details Based on the food safety literature, for given values of \code{c}, \code{r} and \code{t}, the probability of detection in a primary increment is given by, \eqn{p_d=P(X > m)=1-P_{distribution}(X \le m|\mu ,\sigma)} and acceptance probability in \code{t} selected sample is given by \eqn{P_a=P_{binomial}(X \le c|t,p_d)}.
##'
##' If Y be the sum of correlated and identically distributed lognormal random variables X, then the approximate distribution of Y is lognormal
##' distribution with mean \eqn{\mu_y}, standard deviation \eqn{\sigma_y} (see \href{https://doi.org/10.1109/ICC.2006.255040}{Mehta et al (2006)}) where \eqn{E(Y)=mE(X)} and \eqn{V(Y)=mV(X)+cov(X_i,X_j)} for all \eqn{i \ne j =1 \cdots r}.
##'
##' The parameters \eqn{\mu_y} and \eqn{\sigma_y} of the grab sample unit Y is given by,
##' \deqn{\mu_y =\log_{10}{(E[Y])} - {{\sigma_y}^2}/2 \log_e(10) }
##' (see \href{https://doi.org/10.1016/j.foodcont.2013.02.021}{Mussida et al (2013)}). For this package development, we have used fixed \eqn{\sigma_y} value with default value 0.8.
##' @references
##' \itemize{
##' \item Mussida, A., Vose, D. & Butler, F. Efficiency of the sampling plan for {C}ronobacter spp. assuming a Poisson lognormal distribution of the bacteria in powder infant formula and the implications of assuming a fixed within and between-lot variability, Food Control, Elsevier, 2013 , 33 , \href{https://doi.org/10.1016/j.foodcont.2013.02.021}{174-185}.
##' \item Mehta, N.B, Molisch, A.F, Wu, J, &  Zhang, J., 'Approximating the Sum of Correlated Lognormal or, Lognormal-Rice Random Variables,' 2006 IEEE International Conference on Communications, Istanbul, 2006, pp. \href{https://doi.org/10.1109/ICC.2006.255040}{1605-1610}.
##' }
##' @examples
##'   c <-  0
##'   r <-  25
##'   t <-  30
##'   mu <-  -3
##'   distribution <- 'Poisson lognormal'
##'   prob_accept(c, r, t, mu, distribution)
##' @usage prob_accept(c, r, t, mu, distribution, K, m, sd)
##' @export
## we have used the notation as P_D in the paper Quantitative risk assessment for grab sampling inspection of powdered products
prob_accept <- function(c, r, t, mu, distribution, K = 0.25, m = 0, sd = 0.8) {
    p_defect <- NULL
    prob_accept <- NULL
    if (distribution == "Poisson gamma") {
        lambda <- 10^(mu + (sd^2/2) * log(10, exp(1)))
        p_defect <- 1 - extraDistr::pgpois(m, shape = K, rate = K/(lambda * r))
        prob_accept <- stats::pbinom(c, t, p_defect)
    } else if (distribution == "Lognormal") {
        lambda <- 10^(mu + (sd^2/2) * log(10, exp(1)))
        p_defect <- 1 - stats::plnorm(10^m, mu + log(r, 10), sd)
        prob_accept <- stats::pbinom(c, t, p_defect)
    } else if (distribution == "Poisson lognormal") {
        lambda <- 10^(mu + (sd^2/2) * log(10, exp(1)))
        # we borrowed from VGAM package
        dpolono_1 <- function(x, meanlog = 0, sdlog = 1, bigx = 170, ...) {
            mapply(function(x, meanlog, sdlog, ...) {
                if (abs(x) > floor(x)) {
                  0
                } else if (x == Inf) {
                  0
                } else if (x > bigx) {
                  z <- (log(x) - meanlog)/sdlog
                  (1 + (z^2 + log(x) - meanlog - 1)/(2 * x * sdlog^2)) * exp(-0.5 * z^2)/(sqrt(2 * pi) * sdlog * x)
                } else stats::integrate(function(t) exp(t * x - exp(t) - 0.5 * ((t - meanlog)/sdlog)^2), lower = -Inf, upper = Inf, ...)$value/(sqrt(2 *
                  pi) * sdlog * exp(lgamma(x + 1)))
            }, x, meanlog, sdlog, ...)
        }
        ppolono_1 <- function(q, meanlog = 0, sdlog = 1, isOne = 1 - sqrt(.Machine$double.eps), ...) {
            .cumprob <- rep_len(0, length(q))
            .cumprob[q == Inf] <- 1
            q <- floor(q)
            ii <- -1
            while (any(xActive <- ((.cumprob < isOne) & (q > ii)))) .cumprob[xActive] <- .cumprob[xActive] + dpolono_1(ii <- (ii + 1), meanlog, sdlog,
                ...)
            .cumprob
        }

        sam <- function(mu) {
            result <- 1 - ppolono_1(10^m, mu + log(r, 10), sd)
            return(result)
        }
        p_defect <- as.numeric(lapply(mu, sam))
        prob_accept <- stats::pbinom(c, t, p_defect)
    } else {
        warning("please choose one of the given distribution with case sensitive such as 'Poisson' or 'Poisson gamma' or 'Lognormal' or 'Poisson lognormal'")
    }
    return(prob_accept)
}




