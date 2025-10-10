#' Calculate the sampling variance of the log coefficient of variation ratio (lnCVR)
#'
#' This function calculates the sampling variance of the log coefficient of variation ratio (lnCVR),
#' which is used in meta-analysis to compare the variability of a trait between two groups.
#'
#' @param mean.c The mean of the trait in the control group.
#' @param sd.c The standard deviation of the trait in the control group.
#' @param n.c The sample size of the control group.
#' @param cor.c The correlation between the log mean and log standard deviation in the control group.
#' @param mean.e The mean of the trait in the experimental group.
#' @param sd.e The standard deviation of the trait in the experimental group.
#' @param n.e The sample size of the experimental group.
#' @param cor.e The correlation between the log mean and log standard deviation in the experimental group.
#'
#' @return The sampling variance of the lnCVR.
#' @export
#'
#' @examples
#' s2.lnCVR(mean.c = 10, sd.c = 2, n.c = 30, cor.c = 0.8,
#'          mean.e = 12, sd.e = 2.5, n.e = 30, cor.e = 0.85)
s2.lnCVR <- function(mean.c, sd.c, n.c, cor.c, mean.e, sd.e, n.e, cor.e){
  a.c = (sd.c^2)/((n.c)*mean.c^2)
  b.c = 1/(2*(n.c-1))
  c.c = 2*cor.c*(sqrt(a.c*b.c))

  a.e = (sd.e^2)/((n.e)*mean.e^2)
  b.e = 1/(2*(n.e-1))
  c.e = 2*cor.e*(sqrt(a.e*b.e))

  return(a.c + b.c - c.c + a.e + b.e - c.e)
}