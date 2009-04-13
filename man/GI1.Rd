\name{GI1}
\alias{GI1}
\docType{data}
\title{MLGI Output From Study 1}
\description{Output from function MLGI, a matrix of the variance estimates from study 1 of Hamada and Weerandi}
\usage{data(GI1)}
\format{
  A matrix with 17 rows on the following 4 variables.
  \describe{
    \item{\code{ML}}{a numeric vector of maximum likelihood estimates}
    \item{\code{Lower GI}}{a numeric vector of Generalized Inference 95 percent Lower confidence limits}
    \item{\code{Mid GI}}{a numeric vector of Generalized Inference midpoint values}
    \item{\code{Upper GI}}{a numeric vector of Generalized Inference 95 percent Upper confidence limits}
  }
}
\details{GI1 is the output from function MLGI using study 1 data of Hamada and Weerandi. 
It gives the restricted maximum likelihood estimates of variance and the confidence intervals using the method of 
generalized inference as discussed in Hamada and Weerandi. This object can be used to demonstrate the graphing functions.}
\source{Study 1 data is from Michael Hamada and Sam Weerandi (July 2000) Measurement System Assessment 
Via Generalized Inference, Journal of Quality Technology 32, 241-253}
\examples{
data(GI1)
data(sat1)
ratiograph(sat1,GI1)
RRgraph(sat1,GI1)
}
\keyword{datasets}
