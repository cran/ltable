\name{SimData}
\alias{SimData}
\docType{data}
\title{
Simulated interval censored survival data to model hazard.
}
\description{
Simulated interval censored survival data. For details see user guide \emph{ltable2.0.3.pdf}. 
}

\usage{data(SimData)}
\format{
  Description of data frames:
  \describe{
  \var{Fields:}
  \tabular{rlll}{
    1. \tab T \tab Treatment status: 0-1 \tab numeric \cr
    2. \tab C \tab Comorbidity free status: 0-1 \tab numeric \cr
    3. \tab Year \tab Year of event: 1-3 \tab numeric \cr
    4. \tab Counts \tab Generated number of events \tab numeric \cr
    5. \tab Year2 \tab Dummy var. of Year2: 0-1 \tab numeric \cr
    6. \tab Year3 \tab Dummy var. of Year3: 0-1 \tab numeric \cr
    7. \tab offset \tab Person-years of risk \tab numeric \cr
                }
  }
}


\examples{
require(ltable)
data(SimData, package="ltable")
res3<-MCLogLin(Counts~Year2+Year3+T+C, offset=offset, data=SimData, draw=1500, burnin=500)
}
\keyword{datasets}
