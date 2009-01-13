\name{sat2}
\alias{sat2}
\docType{data}
\title{anovasat Output From Study 2}
\description{Gauge R and R output from function anovasat, a matrix of the estimates from study 2 of Hamada and Weerandi}
\usage{data(sat2)}
\format{
  A matrix with 17 rows on the following 3 variables.
  \describe{
    \item{\code{ANOVA est}}{a numeric vector of ANOVA estimates}
    \item{\code{Lower STH}}{a numeric vector of Lower Satterthwaite 95 percent confidence limits}
    \item{\code{Upper STH}}{a numeric vector of Upper Satterthwaite 95 percent confidence limits}
  }
}
\details{sat2 is the Gauge R and R output from function anovasat using study 2 data of Hamada and Weerandi. It gives the ANOVA estimates of variance and the confidence intervals for repeat, reproducibility and gauge using the Satterthwaite approximation as shown on page 244 of Hamada and Weerandi. This can be used to demonstrate the graphing functions.}
\source{Study 1 data is from Michael Hamada and Sam Weerandi (July 2000) Measurement System Assessment Via Generalized Inference, Journal of Quality Technology 32, 241-253}
\examples{
data(sat2)
data(GI2)
ratiograph(sat2,GI2)
RRgraph(sat2,GI2)
}
\keyword{datasets}
