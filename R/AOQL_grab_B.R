##' \code{\link{AOQL_grab_B}} provides the AOQ curve and calculates AOQL value based on average microbial counts.
##' @title Construction of  AOQ curve and calculate AOQL value based on average microbial counts
##' @param c acceptance number
##' @param r number of primary increments in a grab sample or grab sample size
##' @param t number of grab samples
##' @param distribution what suitable microbiological distribution we have used such as  \code{'Poisson gamma'} or \code{'Lognormal'}or \code{'Poisson lognormal'}
##' @param K dispersion parameter of the Poisson gamma distribution (default value 0.25)
##' @param m microbiological limit with default value zero, generally expressed as number of microorganisms in specific sample weight
##' @param sd standard deviation of the lognormal and Poisson-lognormal distributions on the log10 scale (default value 0.8)
##' @param llim the upper limit for graphing the arithmetic mean of cell count
##' @details  Since \eqn{P_a} is the probability of acceptance, \eqn{\lambda} is the arithmetic mean of cell count and the outgoing contaminated arithmetic mean of cell count of primary increments is given by \eqn{AOQ} as the product \eqn{\lambda P_a}.
##'           The quantity \eqn{AOQL} is defined as the maximum proportion of outgoing contaminated primary increments and is given by \deqn{AOQL ={\max_{\lambda \geq 0}}{\lambda P_a}}
##' @seealso  \link{prob_accept}
##' @return AOQ curve and AOQL value based on average microbial counts
##' @examples
##'   c <-  0
##'   r <-  25
##'   t <-  30
##'   distribution <- 'Poisson lognormal'
##'   llim <- 0.20
##'   AOQL_grab_B(c, r, t, distribution, llim)
##' @usage  AOQL_grab_B(c, r, t, distribution,llim, K, m, sd)
##' @export
AOQL_grab_B<- function(c, r, t, distribution,llim, K = 0.25, m = 0, sd = 0.8){
  Sampling_scheme <- NULL  # Initalizing
  P_D <- NULL
  lambda <- seq(0, llim, by = 1e-04)
  mu <- log(lambda, 10)-(sd^2/2)*log(10, exp(1))
  AOQ <- lambda*prob_accept(c, r, t, mu, distribution, K = 0.25, m = 0, sd = 0.8)
  Prob_df <- data.frame(lambda, AOQ)
  f_spr <- function(t, r, c) {
    if (r == 1) {
      sprintf("increments sampling (t=%.0f, c=%.0f)", t, c)
    } else {
      sprintf("grab sampling (t=%.0f, r=%.0f, c=%.0f)", t, r, c)
    }
  }
  Prob <- plyr::rename(Prob_df, c(AOQ = f_spr(t, r, c)))
  melten.Prob <- reshape2::melt(Prob, id = "lambda", variable.name = "Sampling_scheme", value.name = "AOQ")
  ggplot2::ggplot(melten.Prob) + ggplot2::geom_line(ggplot2::aes(x = lambda, y = AOQ, group = Sampling_scheme, colour = Sampling_scheme)) +
    # ggplot2::ggtitle("AOQ curve based on arithmetic mean of cell count") +
    ggplot2::theme_classic() +ggplot2::ylab(expression(AOQ)) + ggthemes::scale_colour_colorblind()+
    ggplot2::xlab(expression("arithmetic mean cell count (" ~ lambda*~")")) +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, size = 10), legend.position = c(0.75, 0.50))+
    # # If we want to add log mean concentration as second x axis
    # ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5), legend.position = c(0.75, 0.50), axis.line.x.top = ggplot2::element_line(color = "red"),
    #                axis.ticks.x.top = ggplot2::element_line(color = "red"), axis.text.x.top = ggplot2::element_text(color = "red"), axis.title.x.top = ggplot2::element_text(color = "red")) +
    # ggplot2::scale_x_continuous(sec.axis = ggplot2::sec_axis(~., name = expression("log mean concentration  (" ~ mu*~")"), breaks = c(0, 0.02, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.50),
    #                                                          labels = c(sprintf("%f", log(0, 10)-(sd^2/2)*log(10, exp(1))),
    #                                                                     sprintf("%f", log(0.02, 10)-(sd^2/2)*log(10, exp(1))),
    #                                                                     sprintf("%f", log(0.05, 10)-(sd^2/2)*log(10, exp(1))),
    #                                                                     sprintf("%f", log(0.10, 10)-(sd^2/2)*log(10, exp(1))),
    #                                                                     sprintf("%f", log(0.15, 10)-(sd^2/2)*log(10, exp(1))),
    #                                                                     sprintf("%f", log(0.20, 10)-(sd^2/2)*log(10, exp(1))),
    #                                                                     sprintf("%f", log(0.25, 10)-(sd^2/2)*log(10, exp(1))),
    #                                                                     sprintf("%f", log(0.30, 10)-(sd^2/2)*log(10, exp(1))),
    #                                                                     sprintf("%f", log(0.50, 10)-(sd^2/2)*log(10, exp(1))))))+
    ggplot2::geom_hline(yintercept=AOQ[which.max(AOQ)],linetype = "dashed")+
    ggplot2::annotate("text", x=4*lambda[which.max(AOQ)], y=AOQ[which.max(AOQ)], label = sprintf("\n AOQL = %0.4f", round(AOQ[which.max(AOQ)], digits = 4)), size=3)
  }






