\name{intsignif}
\alias{intsignif}
\title{P Value of Interaction for Gauge R and R Data}
\description{Displays the P value of the 'part operator' interaction for a gauge R and R data set.}
\usage{intsignif(data)}
\arguments{
  \item{data}{Gauge R and R data set, a matrix with columns named 'resp', 'part' and 'operator'}
}
\details{From a Gauge R and R data set with columns named resp, part and operator, 'intsignif' displays the 
P value of the 'part operator' interaction. The function creates the interaction factor. It dsplays two P values. 
One is the P value of the interaction factor in an ANOVA table. And the other is the likelihood ratio done by doing 
two lmer models, one nested in the other.}
\value{A vector of the two P values}
\references{Package is based on Michael Hamada and Sam Weerandi (July 2000) Measurement System Assessment 
Via Generalized Inference, Journal of Quality Technology 32, 241-253
}
\author{Walter Resch  \email{walt55128@msn.com}}
\note{The argument needs to be a matrix with columns named 'resp', 'part' and 'operator'.}
\seealso{
 \code{ \link{anovasat}},
 \code{ \link{MLGI}},
 \code{ \link{HWstudy1}},
 \code{ \link{HWstudy2}}
}
\examples{
data(HWstudy1)
intsignif(HWstudy1)
}
\keyword{htest}