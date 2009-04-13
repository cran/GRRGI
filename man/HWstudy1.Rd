\name{HWstudy1}
\alias{HWstudy1}
\docType{data}
\title{Gauge R and R: Study 1}
\description{
Gauge R and R study 1 data of Hamada and Weerandi
}
\usage{data(HWstudy1)}
\format{
  A data frame with 120 observations on the following 3 variables.
  \describe{
    \item{\code{part}}{a numeric vector of part numbers}
    \item{\code{operator}}{a numeric vector of operator numbers}
    \item{\code{resp}}{a numeric vector of response data}
  }
}
\details{Gauge R and R data, a matrix of 3 columns and 120 rows. As required for the functions 
anovasat and MLGI, the columns are named: resp, part and operator. There are 20 parts, 3 operators 
and each has two repeats.}
\source{Michael Hamada and Sam Weerandi (July 2000) Measurement System Assessment Via 
Generalized Inference, Journal of Quality Technology 32, 241-253}
\examples{
data(HWstudy1)
anovasat(HWstudy1,0.95)
MLGI(HWstudy1,0.95,1000)
intsignif(HWstudy1)
}
\keyword{datasets}
