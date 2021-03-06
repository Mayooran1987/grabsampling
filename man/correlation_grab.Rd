% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/correlation_grab.R
\name{correlation_grab}
\alias{correlation_grab}
\title{Serial correlation between grab samples}
\usage{
correlation_grab(r, p, d)
}
\arguments{
\item{r}{number of primary increments in a grab sample or grab sample size}

\item{p}{limiting fraction or proportion of contaminated increments}

\item{d}{serial correlation of contamination between the primary increments}
}
\value{
Serial correlation between grab samples
}
\description{
This function calculates the resulting serial correlation between grab samples each having \code{r} primary increments with original serial correlation \code{d}.
}
\details{
The serial correlation between blocks (grab samples) is given by \eqn{d_g} as
\deqn{d_g = [dp(1-p(1-d))^{r-1}]/p_d}
where \eqn{p_d} is the probability of detection in any of the block (grab sample) which is calculated by using \link{prob_detect_single_grab}.
}
\examples{
r <-  25
p <-  0.005
d <-  0.99
correlation_grab(r, p, d)
}
\seealso{
\link{prob_detect_single_grab}
}
