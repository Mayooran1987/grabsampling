##' \code{\link{oc_plan}} provides the Operating Characteristic (OC) curve for known microbiological distribution such as lognormal. The probability of acceptance is plotted against mean log10 concentration.
##' @title Construction of  Operating Characteristic (OC) curve
##' @param c acceptance number
##' @param r number of primary increments in a grab sample or grab sample size
##' @param t number of grab samples
##' @param distribution what suitable distribution we have used such as   \code{'Poisson gamma'} or \code{'Lognormal'} or \code{'Poisson lognormal'}
##' @param sd standard deviation of the lognormal and Poisson-lognormal distributions on the log10 scale (default value 0.8)
##' @param K dispersion parameter of the Poisson gamma distribution (default value 0.25)
##' @param m microbiological limit with default value zero, generally expressed as number of microorganisms in specific sample weight
##' @details Based on the food safety literature, mean concentration is given by \eqn{\lambda = 10^{\mu+log(10)\sigma^2/2}}.
##' @return Operating Characteristic (OC) curve
##' @seealso  \link{prob_accept}
##' @examples
##'   c <-  0
##'   r <-  25
##'   t <-  30
##'   distribution <- 'Poisson lognormal'
##'   oc_plan(c, r, t, distribution)
##' @usage  oc_plan(c, r, t, distribution, K, m, sd)
##' @export
oc_plan <- function(c, r, t, distribution, K = 0.25, m = 0, sd = 0.8) {
    Sampling_scheme <- NULL
    P_a <- NULL
    if (distribution == "Poisson gamma") {
        mu <- seq(-6, -1, 0.01)
        lambda <- 10^(mu + (sd^2/2) * log(10, exp(1)))
        p_a1 <- prob_accept(c, r, t, mu, distribution, K, m, sd)
        Prob_df <- data.frame(mu, p_a1)
        f_spr <- function(t, r, c) {
            if (r == 1) {
                sprintf("primary increment(t=%.0f, c=%.0f)", t, c)
            } else {
                sprintf("grab(t=%.0f, r=%.0f, c=%.0f)", t, r, c)
            }
        }
        Prob <- plyr::rename(Prob_df, c(p_a1 = f_spr(t, r, c)))
        melten.Prob <- reshape2::melt(Prob, id = "mu", variable.name = "Sampling_scheme", value.name = "P_a")
        ggplot2::ggplot(melten.Prob) + ggplot2::geom_line(ggplot2::aes(x = mu, y = P_a, group = Sampling_scheme, colour = Sampling_scheme)) + ggplot2::ggtitle("OC curve based on Poisson gamma distribution") +
            ggplot2::theme_classic() + ggplot2::xlab(expression("log mean concentration  (" ~ mu*~")")) + ggplot2::ylab(expression(P[a])) + ggthemes::scale_colour_colorblind() +
            ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5), legend.position = c(0.85, 0.85), axis.line.x.top = ggplot2::element_line(color = "red"),
                axis.ticks.x.top = ggplot2::element_line(color = "red"), axis.text.x.top = ggplot2::element_text(color = "red"), axis.title.x.top = ggplot2::element_text(color = "red")) +
            ggplot2::scale_x_continuous(sec.axis = ggplot2::sec_axis(~., name = "mean concentration (cfu/g)", breaks = c(-5, -4, -3, -2, -1, 0, 1, 2),
                labels = c(sprintf("%f", 10^(-5 + (sd^2/2) * log(10, exp(1)))), sprintf("%f", 10^(-4 + (sd^2/2) * log(10, exp(1)))), sprintf("%f", 10^(-3 +
                  (sd^2/2) * log(10, exp(1)))), sprintf("%f", 10^(-2 + (sd^2/2) * log(10, exp(1)))), sprintf("%f", 10^(-1 + (sd^2/2) * log(10, exp(1)))),
                  sprintf("%f", 10^(0 + (sd^2/2) * log(10, exp(1)))), sprintf("%f", 10^(1 + (sd^2/2) * log(10, exp(1)))), sprintf("%f", 10^(2 + (sd^2/2) *
                    log(10, exp(1)))))))
    } else if (distribution == "Lognormal") {
        mu <- seq(-3, 0, 0.01)
        lambda <- 10^(mu + (sd^2/2) * log(10, exp(1)))
        p_a1 <- prob_accept(c, r, t, mu, distribution, K, m, sd)
        Prob_df <- data.frame(mu, p_a1)
        f_spr <- function(t, r, c) {
            if (r == 1) {
                sprintf("primary increment(t=%.0f, c=%.0f)", t, c)
            } else {
                sprintf("grab(t=%.0f, r=%.0f, c=%.0f)", t, r, c)
            }
        }
        Prob <- plyr::rename(Prob_df, c(p_a1 = f_spr(t, r, c)))
        melten.Prob <- reshape2::melt(Prob, id = "mu", variable.name = "Sampling_scheme", value.name = "P_a")
        ggplot2::ggplot(melten.Prob) + ggplot2::geom_line(ggplot2::aes(x = mu, y = P_a, group = Sampling_scheme, colour = Sampling_scheme)) + ggplot2::ggtitle("OC curve based on Lognormal distribution") +
            ggplot2::theme_classic() + ggplot2::xlab(expression("log mean concentration  (" ~ mu*~")")) + ggplot2::ylab(expression(P[a])) + ggthemes::scale_colour_colorblind() +
            ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5), legend.position = c(0.85, 0.85), axis.line.x.top = ggplot2::element_line(color = "red"),
                axis.ticks.x.top = ggplot2::element_line(color = "red"), axis.text.x.top = ggplot2::element_text(color = "red"), axis.title.x.top = ggplot2::element_text(color = "red")) +
            ggplot2::scale_x_continuous(sec.axis = ggplot2::sec_axis(~., name = "mean concentration (cfu/g)", breaks = c(-5, -4, -3, -2, -1, 0, 1, 2),
                labels = c(sprintf("%f", 10^(-5 + (sd^2/2) * log(10, exp(1)))), sprintf("%f", 10^(-4 + (sd^2/2) * log(10, exp(1)))), sprintf("%f", 10^(-3 +
                  (sd^2/2) * log(10, exp(1)))), sprintf("%f", 10^(-2 + (sd^2/2) * log(10, exp(1)))), sprintf("%f", 10^(-1 + (sd^2/2) * log(10, exp(1)))),
                  sprintf("%f", 10^(0 + (sd^2/2) * log(10, exp(1)))), sprintf("%f", 10^(1 + (sd^2/2) * log(10, exp(1)))), sprintf("%f", 10^(2 + (sd^2/2) *
                    log(10, exp(1)))))))
    } else if (distribution == "Poisson lognormal") {
        mu <- seq(-6, 0, 0.01)
        lambda <- 10^(mu + (sd^2/2) * log(10, exp(1)))
        p_a1 <- prob_accept(c, r, t, mu, distribution, K, m, sd)
        Prob_df <- data.frame(mu, p_a1)
        f_spr <- function(t, r, c) {
            if (r == 1) {
                sprintf("primary increment(t=%.0f, c=%.0f)", t, c)
            } else {
                sprintf("grab(t=%.0f, r=%.0f, c=%.0f)", t, r, c)
            }
        }
        Prob <- plyr::rename(Prob_df, c(p_a1 = f_spr(t, r, c)))
        melten.Prob <- reshape2::melt(Prob, id = "mu", variable.name = "Sampling_scheme", value.name = "P_a")
        ggplot2::ggplot(melten.Prob) + ggplot2::geom_line(ggplot2::aes(x = mu, y = P_a, group = Sampling_scheme, colour = Sampling_scheme)) + ggplot2::ggtitle("OC curve based on Poisson Lognormal distribution") +
            ggplot2::theme_classic() + ggplot2::xlab(expression("log mean concentration  (" ~ mu*~")")) + ggplot2::ylab(expression(P[a])) + ggthemes::scale_colour_colorblind() +
            ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5), legend.position = c(0.85, 0.85), axis.line.x.top = ggplot2::element_line(color = "red"),
                axis.ticks.x.top = ggplot2::element_line(color = "red"), axis.text.x.top = ggplot2::element_text(color = "red"), axis.title.x.top = ggplot2::element_text(color = "red")) +
            ggplot2::scale_x_continuous(sec.axis = ggplot2::sec_axis(~., name = expression("arithmetic mean cell count (" ~ lambda*~")"), breaks = c(-8, -7, -6, -5, -4, -3, -2,
                -1, 0, 1, 2), labels = c(sprintf("%f", 10^(-8 + (sd^2/2) * log(10, exp(1)))),
                                         sprintf("%f", 10^(-7 + (sd^2/2) * log(10, exp(1)))),
                                         sprintf("%f", 10^(-6 + (sd^2/2) * log(10, exp(1)))),
                                         sprintf("%f", 10^(-5 + (sd^2/2) * log(10, exp(1)))),
                                         sprintf("%f", 10^(-4 + (sd^2/2) * log(10, exp(1)))),
                                         sprintf("%f", 10^(-3 + (sd^2/2) * log(10, exp(1)))),
                                         sprintf("%f", 10^(-2 + (sd^2/2) * log(10, exp(1)))),
                                         sprintf("%f", 10^(-1 + (sd^2/2) * log(10, exp(1)))),
                                         sprintf("%f", 10^(0  + (sd^2/2) * log(10, exp(1)))),
                                         sprintf("%f", 10^(1  + (sd^2/2) * log(10, exp(1)))),
                                         sprintf("%f", 10^(2  + (sd^2/2) * log(10, exp(1)))))))
    } else {
        warning("please choose the one of the given distribution with case sensitive such as 'Poisson gamma' or 'Lognormal' or 'Poisson lognormal'")
    }
}
