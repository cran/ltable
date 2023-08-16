\name{SMdata}
\alias{SMdata}
\docType{data}
\title{
Heart surgery data to model standardized (mortality) ratio.
}
\description{
Counts of patient deaths following heart transplant surgery in 131 hospitals in the US between October 1987 and December 1989. These were analysed by Christiansen and Morris (1996, 1997). Data is also analysed by Peter D. Congdon. Data supplied with file \emph{offsetdata.rda}.
}

\usage{data(SMdata)}
\format{
  Description of data frames:
  \describe{
  \var{Fields:}
  \tabular{rlll}{
    1. \tab y \tab Number of Deaths \tab numeric \cr
    2. \tab o \tab Number of Expected Deaths \tab numeric \cr
                }
          }
      }



\examples{
require(ltable)
data(SMdata, package="ltable")
res2<-MCLogLin(y~1, offset=o, data=SMdata, draw=1500, burnin=500, DIC=TRUE)
}
\references{
Peter D. Congdon. Bayesian Hierarchical Models With Applications Using R (2020) Second Edition. Example 4.5 Hospital Mortality, p.125-26.}

\keyword{datasets}


