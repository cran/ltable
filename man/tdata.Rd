\name{tdata}
\alias{tdata}
\docType{data}
\title{
Tromboembolism Data.
}
\description{
Case-control data first considered by Worcester, J (1971). The data cross-classify tromboembolism and control patients by two risk factors: oral contraceptive user and smoking. Test quantifies boosting effect of contraceptive on odds of tromboembolism. \cr Data used in examples of power analysis.
}

\usage{data(tdata)}
\format{
A grouped data frame with 8 rows of factors' levels combinations. Factors are: smoking status (Yes, No), contraceptive usage (Yes, No), thromboembolism status (Trombol, Control).
  \describe{
    \item{\code{smoker}}{a character vector}
    \item{\code{contraceptive}}{a character vector}
    \item{\code{tromb}}{a character vector}
    \item{\code{Counts}}{a numeric vector}
  }
}
\details{
One can use tables created by function \emph{table_f} transformed with function \emph{tableToData} to appropriate data.frame format with fields of factor, character, logical, and numeric classes. Or one can build data by hand with \emph{data.frame} facility.
}


\examples{
data(tdata, package="ltable")
}

\references{Worcester, J (1971). The relative odds in the 2 by 3 contingency table. American Journal of Epidemiology, 93, 145-149.} 

\keyword{datasets}
