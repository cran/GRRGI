\name{ratiograph}
\alias{ratiograph}
\title{Line Plot of Generalized Inference Confidence Limits for Gauge R and R Ratios}
\description{ratiograph generates a line plot of confidence intervals for gauge R and R ratios from 'MLGI' (generalized inference). It also displays the ANOVA estimate of the ratio from 'anovasat'. }
\usage{ratiograph(sat1,GI1)}
\arguments{
  \item{sat1}{output from 'anovasat'}
  \item{GI1}{output from 'MLGI'}
}
\details{ratiograph generates a line plot of confidence intervals for gauge R and R ratios from 'MLGI' (generalized inference). It also displays the ANOVA estimate of the ratio from 'anovasat'. The data are contained in two matirces that were generated from the functions anovasat and MLGI.}
\value{  Line Plot}
\references{Based on paper by Michael Hamada and Sam Weerandi (July 2000) Measurement System Assessment Via Generalized Inference, Journal of Quality Technology 32, 241-253}
\author{Walter Resch  \email{walt55128@msn.com}}
\note{The arguments need to be the specific outputs from anovasat and MLGI that were run on the same data and at the same percent confidence limits.}
\seealso{
 \code{ \link{anovasat}},
 \code{ \link{MLGI}},
 \code{ \link{sat1}},
 \code{ \link{GI1}},
 \code{ \link{RRgraph}}
}
\examples{
data(sat1)
data(GI1)
ratiograph(sat1,GI1)
}
\keyword{hplot}