\name{sat1}
\alias{sat1}
\docType{data}
\title{anovasat Output From Study 1}
\description{Output from function anovasat, a matrix of the estimates from study 1 of Hamada and Weerandi}
\usage{data(sat1)}
\format{
  A matrix with 17 rows on the following 3 variables.
  \describe{
    \item{\code{ANOVA est}}{a numeric vector of ANOVA estimates}
    \item{\code{Lower STH}}{a numeric vector of Satterthwaite 95 percent Lower confidence limits}
    \item{\code{Upper STH}}{a numeric vector of Satterthwaite 95 percent Upper confidence limits}
  }
}
\details{sat1 is the output from function anovasat using study 1 data of Hamada and Weerandi. 
It gives the ANOVA estimates of variance and the confidence intervals for repeat, reproducibility and 
gauge using the Satterthwaite approximation as presented in Hamada and Weerandi. This object can be used 
to demonstrate the graphing functions.}
\source{Study 1 data is from Michael Hamada and Sam Weerandi (July 2000) Measurement System Assessment Via 
Generalized Inference, Journal of Quality Technology 32, 241-253}
\examples{
data(sat1)
data(GI1)
ratiograph(sat1,GI1)
RRgraph(sat1,GI1)
}
\keyword{datasets}
