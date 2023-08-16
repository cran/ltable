\name{BCdata}
\alias{BCdata}
\docType{data}
\title{
Breast cancer data to model risk
}
\description{
 Breast cancer rates in Iceland by year of birth (11 cohorts from 1840-1849 to 1940-1949) and by age (13 groups from 20-24 to 80-84 years). Analysed by Breslow and Clayton (1993). Data is used also in BUGS Example \dQuote{Ice: non-parametric smoothing in an age-cohort model}. Data supplied with file \emph{offsetdata.rda}.
            }
            
\usage{data(BCdata)}
\format{
   Description of data frames:
  \describe{
  \var{Fields:}
  \tabular{rlll}{
    1. \tab age \tab Age group: 1-13  \tab numeric \cr
    2. \tab year \tab Birth cohort: 1-11 \tab numeric \cr
    3. \tab cases \tab Breast cancer counts \tab numeric \cr
    4. \tab pyr \tab Person-years of risk \tab numeric \cr
                }
          }
      }

\examples{
require(ltable)
data(BCdata, package="ltable")
res1<-MCLogLin(cases~age*year, offset=pyr, data=BCdata, draw=1500, burnin=500)
}

\references{
BUGS. Examples Volume 2. \url{https://www.mrc-bsu.cam.ac.uk/wp-content/uploads/WinBUGS_Vol2.pdf } 
}
\keyword{datasets}
