##' This function produces overlaid Operating Characteristic (OC) curves for any three systematic/random sampling schemes for specified parameters.
##' @title Comparison based on OC curve
##' @param c1,c2,c3 acceptance numbers
##' @param r1,r2,r3 number of primary increments in a grab sample or grab sample size
##' @param t1,t2,t3 number of grab samples
##' @param distribution what distribution we have used such as   \code{'Poisson gamma'} or \code{'Lognormal'}or \code{'Poisson lognormal'}
##' @param K dispersion parameter of the Poisson gamma distribution (default value 0.25)
##' @param m microbiological limit with default value zero, generally expressed as number of microorganisms in specific sample weight
##' @param sd standard deviation of the lognormal and Poisson-lognormal distributions on the log10 scale (default value 0.8)
##' @return overlaid OC curves
##' @seealso \link{prob_accept}
##' @examples
##' c1 <- 0
##' c2 <- 0
##' c3 <- 0
##' r1 <- 25
##' r2 <- 50
##' r3 <- 75
##' t1 <- 10
##' t2 <- 10
##' t3 <- 10
##' distribution <- 'Poisson lognormal'
##' compare_plans_oc(c1, c2, c3, r1, t1, r2, t2, r3, t3, distribution)
##' @usage compare_plans_oc(c1, c2, c3, r1, t1, r2, t2, r3, t3, distribution, K, m, sd)
##' @export
compare_plans_oc <- function(c1, c2, c3, r1, t1, r2, t2, r3, t3, distribution, K = 0.25, m = 0, sd = 0.8) {
    Sampling_scheme <- NULL
    P_a <- NULL
    mu <- seq(-6, 0, 0.01)
    lambda <- 10^(mu + (sd^2/2) * log(10, exp(1)))
    p_a1 <- prob_accept(c1, r1, t1, mu, distribution, K, m, sd)
    p_a2 <- prob_accept(c2, r2, t2, mu, distribution, K, m, sd)
    p_a3 <- prob_accept(c3, r3, t3, mu, distribution, K, m, sd)
    Prob_df <- data.frame(mu, p_a1, p_a2, p_a3)
    f_spr <- function(t, r, c) {
        if (r == 1) {
            sprintf("increments sampling (t=%.0f, c=%.0f)", t, c)
            } else {
                sprintf("grab sampling (t=%.0f, r=%.0f, c=%.0f)", t, r, c)
            }
        }
    Prob <- plyr::rename(Prob_df, c(p_a1 = f_spr(t1, r1, c1), p_a2 = f_spr(t2, r2, c2), p_a3 = f_spr(t3, r3, c3)))
    melten.Prob <- reshape2::melt(Prob, id = "mu", variable.name = "Sampling_scheme", value.name = "P_a")
    ggplot2::ggplot(melten.Prob) + ggplot2::geom_line(ggplot2::aes(x = mu, y = P_a, group = Sampling_scheme, colour = Sampling_scheme)) + ggplot2::xlab(expression("log mean concentration  (" ~ mu*~")")) +
        ggplot2::ylab(expression(P[a])) + ggplot2::theme_classic() + ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, size = 10), legend.position = c(0.75,0.75)) +
        ggthemes::scale_colour_colorblind() + ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5), axis.line.x.top = ggplot2::element_line(color = "red"),axis.ticks.x.top = ggplot2::element_line(color = "red"), axis.text.x.top = ggplot2::element_text(color = "red"), axis.title.x.top = ggplot2::element_text(color = "red")) +
        ggplot2::scale_x_continuous(sec.axis = ggplot2::sec_axis(~., name = expression("arithmetic mean cell count (" ~ lambda*~")"), breaks = c(-8, -7, -6, -5, -4, -3, -2, -1, 0, 1, 2), labels = c(sprintf("%f", 10^(-8 + (sd^2/2) * log(10, exp(1)))), sprintf("%f", 10^(-7 + (sd^2/2) * log(10, exp(1)))), sprintf("%f", 10^(-6 + (sd^2/2) * log(10, exp(1)))), sprintf("%f", 10^(-5 + (sd^2/2) * log(10, exp(1)))), sprintf("%f", 10^(-4 + (sd^2/2) * log(10, exp(1)))),
                                                                                                                                                                          sprintf("%f", 10^(-3 + (sd^2/2) * log(10, exp(1)))), sprintf("%f", 10^(-2 + (sd^2/2) * log(10, exp(1)))), sprintf("%f", 10^(-1 + (sd^2/2) * log(10, exp(1)))), sprintf("%f", 10^(0 + (sd^2/2) * log(10, exp(1)))), sprintf("%f", 10^(1 + (sd^2/2) * log(10, exp(1)))), sprintf("%f", 10^(2 + (sd^2/2) * log(10, exp(1)))))))
    }

