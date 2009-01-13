\name{sat1}
\alias{sat1}
\docType{data}
\title{anovasat Output From Study 1}
\description{Gauge R and R output from function anovasat, a matrix of the estimates from study 1 of Hamada and Weerandi}
\usage{data(sat1)}
\format{
  A matrix with 17 rows on the following 3 variables.
  \describe{
    \item{\code{ANOVA est}}{a numeric vector of ANOVA estimates}
    \item{\code{Lower STH}}{a numeric vector of Lower Satterthwaite 95 percent confidence limits}
    \item{\code{Upper STH}}{a numeric vector of Upper Satterthwaite 95 percent confidence limits}
  }
}
\details{sat1 is the Gauge R and R output from function anovasat using study 1 data of Hamada and Weerandi. It gives the ANOVA estimates of variance and the confidence intervals for repeat, reproducibility and gauge using the Satterthwaite approximation as shown on page 244 of Hamada and Weerandi. This can be used to demonstrate the graphing functions.}
\source{Study 1 data is from Michael Hamada and Sam Weerandi (July 2000) Measurement System Assessment Via Generalized Inference, Journal of Quality Technology 32, 241-253}
\examples{
data(sat1)
data(GI1)
ratiograph(sat1,GI1)
RRgraph(sat1,GI1)
}
\keyword{datasets}
