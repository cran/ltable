\name{ltable-package}
\alias{ltable}
\docType{package}
\title{
ltable}

\description{
\bold{Constructs tables of counts and proportions out of data sets.}
}


\details{

Package: ltable

Type: Package

Title: Easy to make (Lazy) tables

Version: 1.0

Date: 2018-08-10

Depends: R (>= 3.4.0)

Suggests: utils, stats

License: GPL (>= 2)

}


\author{
{Ocheredko Oleksandr}
\email{Ocheredko@yahoo.com}
}

\note{
You can:
\enumerate{
\item{construct tables with data set fields of factor, character, logical, and numeric classes;}
\item{insert tables into Excell and Word documents.}
		}
	}


\seealso{
\code{\link[etable:etable]{etable}}
\code{\link[xtable:xtable]{xtable}}
\code{\link[tables:tabular]{tabular}}


}

\examples{
data(sdata, package="ltable")
table_f(sdata, "a")
table_f(sdata, "a", MV=TRUE, extended=TRUE)
table_f(sdata, "a,b,c")
table_f(sdata, "a,b,c,d", type=2, digits=3)
table_f(sdata, "b,c,a,d", MV=TRUE, extended=TRUE, cb=TRUE)
}

\keyword{utilities}
