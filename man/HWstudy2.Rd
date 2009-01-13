\name{HWstudy2}
\alias{HWstudy2}
\docType{data}
\title{Gauge R and R: Study 2}
\description{Gauge R and R study 2 data of Hamada and Weerandi}
\usage{data(HWstudy2)}
\format{
  A data frame with 40 observations on the following 3 variables.
  \describe{
    \item{\code{part}}{a numeric vector of part numbers}
    \item{\code{operator}}{a numeric vector of operator numbers}
    \item{\code{resp}}{a numeric vector of responses}
  }
}
\details{Gauge R and R data, a matrix of 3 columns and 40 rows. As required for the functions anovasat and MLGI, the columns are named: resp, part and operator. There are 10 parts, 2 operators and each has two repeats.}
\source{Michael Hamada and Sam Weerandi (July 2000) Measurement System Assessment Via Generalized Inference, Journal of Quality Technology 32, 241-253}
\examples{
data(HWstudy2)
anovasat(HWstudy2,0.95)
MLGI(HWstudy2,0.95,1000)
intsignif(HWstudy2)
}
\keyword{datasets}
