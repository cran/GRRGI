\name{RRgraph}
\alias{RRgraph}
\title{Comparison Plot of Satterthwaite and Generalized Inference Confidence Limits}
\description{RRgraph generates a line plot comparing the gauge, reproducibility and repeatability confidence intervals from 'anovasat' (Satterthwaite) to the confidence intervals from 'MLGI' (generalized inference).}
\usage{RRgraph(sat1,GI1)}
\arguments{
  \item{sat1}{output from 'anovasat'}
  \item{GI1}{output from 'MLGI'}
}
\details{RRgraph generates a line plot comparing the gauge, reproducibility and repeatability confidence intervals from 'anovasat' (Satterthwaite) to the confidence intervals from 'MLGI' (generalized inference). The confidence intervals are contained in two matirces that were generated from the functions anovasat and MLGI.}
\value{  Line plot}
\references{Based on paper by Michael Hamada and Sam Weerandi (July 2000) Measurement System Assessment Via Generalized Inference, Journal of Quality Technology 32, 241-253}
\author{Walter Resch  \email{walt55128@msn.com}}
\note{The arguments need to be the specific outputs from anovasat and MLGI that were run on the same data and at the same percent confidence limits.}
\seealso{
 \code{ \link{anovasat}},
 \code{ \link{MLGI}},
 \code{ \link{sat1}},
 \code{ \link{GI1}},
 \code{ \link{ratiograph}}
}
\examples{
data(sat1)
data(GI1)
RRgraph(sat1,GI1)
}
\keyword{hplot}
