##' This function allows comparison of different sampling schemes, which can be systematic and random sampling of primary increments or grab sampling of blocks of primary increments.  A graphical display of the probability of detection \eqn{P_D} or probability of non detection \eqn{P_{ND}} versus fraction nonconforming \eqn{p} for up to four selected schemes will be produced.
##' @title Probability of detection or non detection  versus fraction nonconforming curve
##' @param c1,c2,c3,c4 acceptance numbers
##' @param r1,r2,r3,r4 number of primary increments in a grab sample or grab sample size
##' @param t1,t2,t3,t4 number of grab samples
##' @param d serial correlation of contamination between the primary increments
##' @param N length of the production
##' @param method1,method2,method3,method4 what sampling method we have applied such as \code{'systematic'} or \code{'random'} selection methods
##' @param plim the upper limit for graphing the fraction nonconforming or proportion of contaminated increments
##' @param type what type of graph we want to produce such as \code{D} or \code{ND}. \code{\link{compare_plans}} produces a graphical display of \eqn{P_D} or \eqn{P_{ND}} versus \eqn{p} depending on the \code{D} or \code{ND} of type
##' @return Probability of detection or non detection vs limiting fraction curves
##' @examples
##' c1 <- 0
##' c2 <- 0
##' c3 <- 0
##' c4 <- 0
##' r1 <- 1
##' r2 <- 10
##' r3 <- 30
##' r4 <- 75
##' t1 <- 750
##' t2 <- 75
##' t3 <- 25
##' t4 <- 10
##' d <- 0.99
##' N <- 1e9
##' method1 <- method2 <- method3 <- method4 <- 'systematic'
##' plim <- 0.10
##' compare_plans(d, N, plim, type ='D', c1, r1, t1, method1, c2, r2, t2, method2)
##' compare_plans(d, N, plim, type ='D', c1, r1, t1, method1, c2, r2, t2, method2,
##'                         c3, r3, t3, method3)
##' compare_plans(d, N, plim, type ='D', c1, r1, t1, method1, c2, r2, t2, method2,
##'                         c3, r3, t3, method3, c4, r4, t4, method4)
##' compare_plans(d, N, plim, type ='ND', c1, r1, t1, method1, c2, r2, t2, method2,
##'                         c3, r3, t3, method3, c4, r4, t4, method4)
##'
##' @usage compare_plans(d, N, plim, type, c1, r1, t1, method1, c2, r2, t2, method2,
##'                      c3, r3, t3, method3, c4, r4, t4, method4)
##' @export
compare_plans <- function(d, N, plim, type, c1, r1, t1, method1, c2 = NULL, r2 = NULL, t2 = NULL, method2 = NULL, c3 = NULL, r3 = NULL, t3 = NULL, method3 = NULL, c4 = NULL, r4 = NULL, t4 = NULL, method4 = NULL) {
    Sampling_scheme <- NULL  # Initalizing
    P_D <- NULL
    p <- seq(1e-05, plim, by = 1e-05)
    f_spr <- function(t, r, c, method) {
        if (method == "systematic") {
            if (r == 1) {
                sprintf("systematic increments sampling (t=%.0f, c=%.0f)", t, c)
                } else {
                    sprintf("systematic grab sampling (t=%.0f, r=%.0f, c=%.0f)", t, r, c)
                    }
            } else {
                if (r == 1) {
                    sprintf("random increments sampling (t=%.0f, c=%.0f)", t, c)
                    } else {
                        sprintf("random grab sampling (t=%.0f, r=%.0f, c=%.0f)", t, r, c)
                    }
            }
        }
    if(is.null(c4) && is.null(r4) && is.null(t4) ) {
        if(is.null(c3) && is.null(r3) && is.null(t3) ) {
            p_d1 <- prob_detect(c1, r1, t1, d, p, N, method1)
            p_d2 <- prob_detect(c2, r2, t2, d, p, N, method2)
            Prob_df <- data.frame(p, p_d1, p_d2)
            Prob <- plyr::rename(Prob_df, c(p_d1 = f_spr(t1, r1, c1, method1), p_d2 = f_spr(t2, r2, c2,method2)))
        } else {
            p_d1 <- prob_detect(c1, r1, t1, d, p, N, method1)
            p_d2 <- prob_detect(c2, r2, t2, d, p, N, method2)
            p_d3 <- prob_detect(c3, r3, t3, d, p, N, method3)
            Prob_df <- data.frame(p, p_d1, p_d2, p_d3)
            Prob <- plyr::rename(Prob_df, c(p_d1 = f_spr(t1, r1, c1, method1), p_d2 = f_spr(t2, r2, c2, method2), p_d3 = f_spr(t3, r3, c3, method3)))
        }
    } else {
        p_d1 <- prob_detect(c1, r1, t1, d, p, N, method1)
        p_d2 <- prob_detect(c2, r2, t2, d, p, N, method2)
        p_d3 <- prob_detect(c3, r3, t3, d, p, N, method3)
        p_d4 <- prob_detect(c4, r4, t4, d, p, N, method4)
        Prob_df <- data.frame(p, p_d1, p_d2, p_d3, p_d4)
        Prob <- plyr::rename(Prob_df, c(p_d1 = f_spr(t1, r1, c1, method1), p_d2 = f_spr(t2, r2, c2, method2), p_d3 = f_spr(t3, r3, c3, method3), p_d4 = f_spr(t4, r4, c4, method4)))
    }
    melten.Prob <- reshape2::melt(Prob, id = "p", variable.name = "Sampling_scheme", value.name = "P_D")
    if (type == "D") {
        ggplot2::ggplot(melten.Prob) + ggplot2::geom_line(ggplot2::aes(x = p, y = P_D, group = Sampling_scheme, colour = Sampling_scheme)) + ggplot2::ylab(expression(P[D])) +
            ggplot2::theme_classic() + ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, size = 10), legend.position = c(0.75, 0.25)) + ggthemes::scale_colour_colorblind()
        } else if (type == "ND"){
            ggplot2::ggplot(melten.Prob) + ggplot2::geom_line(ggplot2::aes(x = p, y = 1-P_D, group = Sampling_scheme, colour = Sampling_scheme)) + ggplot2::ylab(expression(P[ND])) +
                ggplot2::theme_classic() + ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, size = 10), legend.position = c(0.75, 0.75)) + ggthemes::scale_colour_colorblind()
        }
    }
