\name{GRRGI-package}
\alias{GRRGI-package}
\alias{GRRGI}
\docType{package}
\title{Components of Variance and their Confidence Intervals for Gauge R and R Studies
}
\description{Calculates components of Variance for Gauge R and R Studies using ANOVA and Maximum Likelihood. 
Also generates the confidence intervals using the Satterthwaite approximation with the ANOVA, and using the 
method of Generalized Inference with the maximum likelihood. Also generates line plots to compare the Satterthwaite 
and Generalized Inference confidence limits.}
\details{
\tabular{ll}{
Package: \tab GRRGI\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2008-12-28\cr
License: \tab GPL-2\cr
LazyLoad: \tab yes\cr
}

Gauge R and R data need to be in a matrix with columns named 'resp', 'part' and 'operator'. Example of data sets are provided in 'HWstudy1' and 'HWstudy2'. 

'anovasat' calculates the components of variance using ANOVA and the the confidence intervals for repeatability, reproducibility and total gauge using the Satterthwaite approximation.

'MLGI' calculates the components of Variance using maximum likelihood and the the confidence intervals using the method of generalized inference.
}
\author{Walter Resch  \email{walt55128@msn.com}}
\references{Package is based on Michael Hamada and Sam Weerandi (July 2000) Measurement System Assessment Via Generalized Inference, Journal of Quality Technology 32, 241-253}
\keyword{ package }
\seealso{
 \code{ \link{MLGI}},
 \code{ \link{anovasat}},
 \code{ \link{intsignif}},
 \code{ \link{RRgraph}},
 \code{ \link{ratiograph}}
}
\examples{
data(HWstudy1)
sat<-anovasat(HWstudy1,0.95)
ML<-MLGI(HWstudy1,0.95,1000)
intsignif(HWstudy1)
RRgraph(sat,ML)
ratiograph(sat,ML)
}
