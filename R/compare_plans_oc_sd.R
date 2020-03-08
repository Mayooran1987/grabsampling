##' This function produces overlaid Operating Characteristic (OC) curves for any three systematic/random sampling schemes for specified parameters with different standard deviation vlues.
##' @title Comparison based on OC curve with different standard deviations
##' @param c1,c2,c3 acceptance numbers
##' @param r1,r2,r3 number of primary increments in a grab sample or grab sample size
##' @param t1,t2,t3 number of grab samples
##' @param distribution what distribution we have used such as   \code{'Poisson gamma'} or \code{'Lognormal'}or \code{'Poisson lognormal'}
##' @param K dispersion parameter of the Poisson gamma distribution (default value 0.25)
##' @param m microbiological limit with default value zero, generally expressed as number of microorganisms in specific sample weight
##' @param sd1,sd2,sd3 standard deviations of the lognormal and Poisson-lognormal distributions on the log10 scale
##' @return overlaid OC curves
##' @seealso \link{prob_accept}
##' @examples
##' c1 <- 0
##' c2 <- 0
##' c3 <- 0
##' r1 <- 25
##' r2 <- 25
##' r3 <- 25
##' t1 <- 30
##' t2 <- 30
##' t3 <- 30
##' sd1 <- 0.2
##' sd2 <- 0.4
##' sd3 <- 0.8
##' distribution <- 'Poisson lognormal'
##' compare_plans_oc_sd(c1, c2, c3, r1, t1, r2, t2, r3, t3, sd1, sd2, sd3, distribution)
##' @usage compare_plans_oc_sd(c1, c2, c3, r1, t1, r2, t2, r3, t3, sd1, sd2, sd3, distribution, K, m)
##' @export
compare_plans_oc_sd <- function(c1, c2, c3, r1, t1, r2, t2, r3, t3, sd1, sd2, sd3, distribution, K = 0.25, m = 0) {
  Sampling_scheme <- NULL
  P_a <- NULL
  mu <- seq(-6, 0, 0.05)
  # lambda <- 10^(mu + (sd^2/2) * log(10, exp(1)))
  p_a1 <- grabsampling::prob_accept(c1, r1, t1, mu, distribution, sd=sd1, K, m)
  p_a2 <- grabsampling::prob_accept(c2, r2, t2, mu, distribution, sd=sd2, K, m)
  p_a3 <- grabsampling::prob_accept(c3, r3, t3, mu, distribution, sd=sd3, K, m)
  Prob_df <- data.frame(mu, p_a1, p_a2, p_a3)
  f_spr <- function(t, r, c, sd) {
    if (r == 1) {
      sprintf("increments sampling (t=%.0f, c=%.0f, sd=%.2f)", t, c, sd)
    } else {
      sprintf("grab sampling (t=%.0f, r=%.0f, c=%.0f, sd=%.2f)", t, r, c, sd)
    }
  }
  Prob <- plyr::rename(Prob_df, c(p_a1 = f_spr(t1, r1, c1, sd1), p_a2 = f_spr(t2, r2, c2, sd2), p_a3 = f_spr(t3, r3, c3, sd3)))
  melten.Prob <- reshape2::melt(Prob, id = "mu", variable.name = "Sampling_scheme", value.name = "P_a")
  ggplot2::ggplot(melten.Prob) + ggplot2::geom_line(ggplot2::aes(x = mu, y = P_a, group = Sampling_scheme, colour = Sampling_scheme)) +
    ggplot2::xlab(expression("log mean concentration  (" ~ mu*~")")) +
    ggplot2::ylab(expression(P[a])) + ggplot2::theme_classic() +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, size = 10),legend.position = c(0.75,0.75)) +
    ggthemes::scale_colour_colorblind()
}
