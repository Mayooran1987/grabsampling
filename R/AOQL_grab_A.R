##' \code{\link{AOQL_grab_A}} provides the AOQ curve and calculates AOQL value based on limiting fraction of contaminated increments.
##' @title Construction of  AOQ curve and calculate AOQL value based on limiting fraction
##' @param c acceptance number
##' @param r nurber of primary increments in a grab sample or grab sample size
##' @param t number of grab samples
##' @param d serial correlation of contamination between the primary increments
##' @param N length of the production
##' @param method what sampling method we have applied such as \code{'systematic'} or \code{'random'} selection methods
##' @param plim the upper limit for graphing the fraction nonconforming or proportion of contaminated increments
##' @details  Since \eqn{P_{ND}} is the probability of non-detection, \eqn{p} is the limiting fraction of contaminated increments and the outgoing contaminated proportion of primary increments is given by \eqn{AOQ} as the product \eqn{pP_{ND}}.
##'           The quantity \eqn{AOQL} is defined as the maximum proportion of outgoing contaminated primary increments and is given by \deqn{AOQL ={\max_{0\leq p\leq 1}}{pP_{ND}}}
##' @seealso  \link{prob_detect}
##' @return AOQ curve and AOQL value based on on limiting fraction
##' @examples
##'   c <-  0
##'   r <-  25
##'   t <-  30
##'   d <-  0.99
##'   N <-  1e9
##'   method <- 'systematic'
##'   plim <- 0.30
##'   AOQL_grab_A(c, r, t, d, N, method, plim)
##' @usage  AOQL_grab_A(c, r, t, d, N, method, plim)
##' @export
AOQL_grab_A <- function(c, r, t, d, N, method, plim){
  Sampling_scheme <- NULL  # Initalizing
  P_D <- NULL
  p <- seq(1e-05, plim, by = 1e-05)
  if (method == "systematic") {
    f_spr <- function(t, r, c) {
      if (r == 1) {
        sprintf("systematic increments sampling (t=%.0f, c=%.0f)", t, c)
      } else {
        sprintf("systematic grab sampling (t=%.0f, r=%.0f, c=%.0f)", t, r, c)
      }
    }
  } else {
    f_spr <- function(t, r, c) {
      if (r == 1) {
        sprintf("random increments sampling (t=%.0f, c=%.0f)", t, c)
      } else {
        sprintf("random grab sampling (t=%.0f, r=%.0f, c=%.0f)", t, r, c)
      }
    }
  }
  AOQ <- p*(1-prob_detect(c, r, t, d, p, N, method))
  Prob_df <- data.frame(p, AOQ)
  Prob <- plyr::rename(Prob_df, c(AOQ = f_spr(t, r, c)))
  melten.Prob <- reshape2::melt(Prob, id = "p", variable.name = "Sampling_scheme", value.name = "AOQ")
  ggplot2::ggplot(melten.Prob) + ggplot2::geom_line(ggplot2::aes(x = p, y = AOQ , group = Sampling_scheme, colour = Sampling_scheme)) +
    # ggplot2::ggtitle("AOQ curve based on limiting fraction of contaminated increments") +
    ggplot2::ylab(expression(AOQ)) +ggplot2::xlab(expression("limiting fraction (" ~ p*~")"))+
    ggplot2::theme_classic() + ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, size = 10), legend.position = c(0.75, 0.50)) + ggthemes::scale_colour_colorblind() +
    ggplot2::geom_hline(yintercept=AOQ[which.max(AOQ)],linetype = "dashed")+
    ggplot2::annotate("text", x=4*p[which.max(AOQ)], y=AOQ[which.max(AOQ)], label = sprintf("\n AOQL = %0.4f", round(AOQ[which.max(AOQ)], digits = 4)), size=3)
}
