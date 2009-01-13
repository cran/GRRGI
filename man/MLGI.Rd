\name{MLGI}
\alias{MLGI}
\title{Maximum Likelihood Estimates of Variance and Generalized Inference Confidence Limits for Gauge R and R Data}
\description{Returns a matrix of the Maximum Likelihood Estimates of Variance for gauge R and R studies. Also provides the confidence intervals using Gerneralized Inference as shown on pages 245-246 in the paper by Hamada and Weerandi.}
\usage{MLGI(data, prob, N)}
\arguments{
  \item{data}{Gauge R and R data, a matrix with columns named 'resp', 'part' and 'operator' in any order}
  \item{prob}{desired confidence limits such as 0.95 - a constant between 0 and 1}
  \item{N}{desired number of simulations for Gerneralized Inference}
}
\details{Returns a matrix with 17 rows and four columns with the estimates: The standard deviations:  'SD:Part', 'SD:Operator', 'SD:PartOp', 'SD:Repeat', 'SD:Reproduce', 'SD:Gauge', 'SD:Total'.  The Variances:  'Var:Part', 'Var:Operator', 'Var:PartOp', 'Var:Repeat', 'Var:Reproduce', 'Var:Gauge', 'Var:Total'.  The ratios:  'Gauge/Total', 'Gauge/Parts', 'Repeat/Gauge'.  The columns are 'ML est', 'Lower GI', 'Mid GI', 'Upper GI'. The Lower and Upper GI are the Generalized Inference confidence limits. 'MLGI' will eliminate the 'part*operator' interaction factor if the P value of the likelihood ratio test is greater than 0.25.}
\value{A matrix with 17 rows and four columns with the estimates.}
\references{Package is based on Michael Hamada and Sam Weerandi (July 2000) Measurement System Assessment Via Generalized Inference, Journal of Quality Technology 32, 241-253}
\author{Walter Resch  \email{walt55128@msn.com}}
\note{The 'data' argument needs to be a matrix with columns named 'resp', 'part' and 'operator'.}
\seealso{
  \code{\link{anovasat}},
  \code{\link{RRgraph}},
  \code{\link{ratiograph}}
}
\examples{
data(HWstudy1)
MLGI(HWstudy1,0.95,1000)
}
\keyword{htest}
