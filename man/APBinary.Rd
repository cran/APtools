\name{APBinary}
\alias{APBinary}
\title{Estimating the AP and the AUC for Binary Outcome Data.}
\description{
  This function calculates the estimates of the AP and AUC
  for binary outcomes as well as their confidence intervals
  using the perturbation or the nonparametric bootstrap
  resampling method.
}
\usage{
  APBinary(status, marker, cut.values = NULL,
    method = "bootstrap", alpha = 0.95, B = 1000)
}
\arguments{
  \item{status}{Binary indicator, 1 indicates case / the
  class of prediction interest and 0 otherwise.}

  \item{marker}{Numeric risk score or marker value. Data
  can be continuous or ordinal.}

  \item{cut.values}{marker values to use as a cut-off for
  calculation of positive predictive values (PPV) and true
  positive fractions (TPF). The default value is NULL.}

  \item{method}{Method to obtain confidence intervals. The
  default value is "bootstrap". Another option is
  "perturbation".}

  \item{alpha}{Confidence level.  The default level is
  0.95.}

  \item{B}{Number of resampling to obtain confidence
  interval. The default value is 1000.}
}
\value{
  an object of class "APBinary" which is a list with
  components:

  \item{ap_summary}{Summary of the AP, including the proportion of
  cases, a point estimate of AP, and their corresponding
  confidence intervals.}

  \item{auc_summary}{Summary of the AUC, including a point
  estimate of AUC with a confidence interval.}

  \item{PPV}{Available object, positive predictive values at the
  unique risk score in the data.}

  \item{TPF}{Available object, true positive fractions at the
  unique risk score in the data.}
}
\references{
  Yuan, Y., Su, W., and Zhu, M. (2015). Threshold-free
  measures for assessing the performance of medical
  screening tests. Frontiers in Public Health, 3.57.

  Bingying Li (2015) Threshold-free Measure for Assessing
  the Performance of Risk Prediction with Censored Data,
  MSc. thesis, Simon Fraser University, Canada
}
\examples{
status=c(rep(1,10),rep(0,1),rep(1,18),rep(0,11),rep(1,25),
	rep(0,44),rep(1,85),rep(0,176))
marker=c(rep(7,11),rep(6,29),rep(5,69),rep(4,261))
cut.values=sort(unique(marker)[-1])
out1 <- APBinary(status,marker,cut.values)
out1$PPV
out1$TPF
out2 <- APBinary(status,marker,method="perturbation",
	alpha=0.90,B=1500)
out2$ap_summary
out2$auc_summary
}

