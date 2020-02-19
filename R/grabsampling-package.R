##' This package provides the probability of detection calculation for grab samples selection by various method of samplings such as systematic or random and also it is useful to generate the comparison curves.
##' Moreover, this package calculates the probability of acceptance calculations based on suitable microbiological distributions such as Poisson gamma or Lognormal or Poisson lognormal and also provides a comparison based on OC curves with different sampling schemes. Most of the researchers have studied the uncorrelated case, but in this study, we have a high spatial correlation between contamination of primary increments. For this package development, we used default standard deviation as 0.8 and also spatial correlation not affected for composite mean for the probability of acceptance calculation. A future version will be included deeply study about variability effects in grab sample selection.
##' @name grabsampling-package
##' @docType package
##' @title Probability of detection for grab sample selection
##' @author Mayooran Thevaraja, Kondaswamy Govindaraju, Mark Bebbington
##' @references
##' \itemize{
##' \item Bhat, U., & Lal, R. (1988). Number of successes in Markov trials. Advances in Applied Probability, 20(3), \href{https://doi.org/10.2307/1427041}{677-680}.
##' \item Jongenburger, I., Besten, H.M., & Zwietering, M.H. (2015). Statistical aspects of food safety sampling. Annual review of food science and technology, 6, \href{https://doi.org/10.1146/annurev-food-022814-015546}{479-503}.
##' \item Mussida, A., Vose, D. & Butler, F. Efficiency of the sampling plan for {C}ronobacter spp. assuming a Poisson lognormal distribution of the bacteria in powder infant formula and the implications of assuming a fixed within and between-lot variability, Food Control, Elsevier, 2013 , 33 , \href{https://doi.org/10.1016/j.foodcont.2013.02.021}{174-185}.
##' \item Van Schothorst, M., Zwietering, M., Ross, T., Buchanan, R. & Cole, M., Relating microbiological criteria to food safety objectives and performance objectives  Food Control , 2009 , 20 , \href{https://doi.org/10.1016/j.foodcont.2008.11.005}{967-979}.
##' }
NULL
