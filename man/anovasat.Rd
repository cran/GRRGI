\name{anovasat}
\alias{anovasat}
\title{ANOVA Estimates of Variance and Satterthwaite Confidence Limits for Gauge R and R Data}
\description{Returns a matrix of the ANOVA Estimates of Variance for gauge R and R studies. Also 
provides the confidence intervals for repeatability, reproducibility and total gauge using the 
Satterthwaite approximation as shown on page 244 of Hamada and Weerandi. }
\usage{anovasat(data, prob)}
\arguments{
  \item{data}{Gauge R and R data, a matrix with columns named 'resp', 'part' and 'operator' in any order}
  \item{prob}{desired confidence level such as 0.95 - a constant between 0 and 1}
}
\details{Returns a matrix with 17 rows and three columns with the estimates:  
SD:Part, SD:Operator, SD:PartOp, SD:Repeat, SD:Reproduce, SD:Gauge, SD:Total (standard deviations), then   
Var:Part, Var:Operator, Var:PartOp, Var:Repeat, Var:Reproduce, Var:Gauge, Var:Total (variances), then  
Gauge/Total, Gauge/Parts, Repeat/Gauge (ratios.)  

The columns are 'ANOVA est' (the ANOVA based variance estimate,) 'Lower STH', 'Upper STH'. 
The 'Lower STH' and 'Upper STH' are the Satterthwaite based confidence limits. 
They are calculated for repeatability, reproducibility and total gauge.  'anovasat' will 
eliminate the 'part operator' interaction factor if its P value is greater than 0.25.  }
\value{  A matrix with 17 rows and three columns with the estimates.}
\references{Package is based on Michael Hamada and Sam Weerandi (July 2000) Measurement 
System Assessment Via Generalized Inference, Journal of Quality Technology 32, 241-253}
\author{Walter Resch  \email{walt55128@msn.com}}
\note{The 'data' argument needs to be a matrix with columns named 'resp', 'part' and 'operator'.}
\seealso{
 \code{ \link{MLGI}},
 \code{ \link{RRgraph}},
 \code{ \link{ratiograph}}
}
\examples{
data(HWstudy1)
anovasat(HWstudy1,0.95)
}
\keyword{htest}
