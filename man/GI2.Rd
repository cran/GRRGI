\name{GI2}
\alias{GI2}
\docType{data}
\title{MLGI Output From Study 2}
\description{Gauge R and R output from function MLGI, a matrix of the variance estimates from study 2 of Hamada and Weerandi}
\usage{data(GI2)}
\format{
  A matrix with 17 rows on the following 4 variables.
  \describe{
    \item{\code{ML}}{a numeric vector of Maximum Likelihood estimates}
    \item{\code{Lower GI}}{a numeric vector of Lower Generalized Inference 95 percent confidence limits}
    \item{\code{Mid GI}}{a numeric vector of midpoint Generalized Inference values}
    \item{\code{Upper GI}}{a numeric vector of Upper Generalized Inference 95 percent confidence limits}
  }
}
\details{GI2 is the Gauge R and R output from function MLGI using study 2 data of Hamada and Weerandi. 
It gives the maximum likelihood estimates of variance and the confidence intervals using the method of 
generalized inference as discussed in Hamada and Weerandi. This can be used to demonstrate the graphing functions.}
\source{Study 2 data is from Michael Hamada and Sam Weerandi (July 2000) Measurement System Assessment Via 
Generalized Inference, Journal of Quality Technology 32, 241-253}
\examples{
data(GI2)
data(sat2)
ratiograph(sat2,GI2)
RRgraph(sat2,GI2)
}
\keyword{datasets}