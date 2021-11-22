\name{ltable-package}
\alias{ltable}
\docType{package}
\title{
ltable}

\description{
\bold{Constructs tables of counts and proportions out of data sets.}
}


\details{

In order to perform log-linear and power analyses for tabulated data GSL: GNU Scientific Library has to be installed first (reference 1). Log-linear analysis features some advantages against glm \{stats\}, first of all due to stability of GSL IWLS algorithms that insures distinctly less biased covariances estimates, pivot issue for implemented power analysis. In some instances hypothesis testing of higher order effects disagrees with that of glm on account of larger GSL estimated errors. Another though related enhancement is distinct better fit assessed by sum of squared differences between observed and expected counts. Results of power analysis backed up with MCMC delivered approach (reference 2).  

}


\author{
{Ocheredko Oleksandr}
\email{Ocheredko@yahoo.com}
}

\note{
You can:
\enumerate{
\item{construct tables with data set fields of factor, character, logical, and numeric classes;}
\item{insert tables into Excel and Word documents using clipboard, into LaTeX, HTML, Markdown and reStructuredText documents by the knitr::kable agency.}
		}
	}


\examples{
require(ltable)
data(sdata, package="ltable")
table_f(sdata, "a")
table_f(sdata, "a", MV=TRUE, extended=TRUE)
table_f(sdata, "a,b,c")
knitr::kable(table_f(sdata, "a,b,c,d", type=2, digits=3))
table_f(sdata, "b,c,a,d", MV=TRUE, extended=TRUE, cb=TRUE)
}
\references{
1. GSL: GNU Scientific Library: \url{http://www.gnu.org/software/gsl/},\cr 2. Ocheredko O.M. MCMC Bootstrap Based Approach to Power and Sample Size Evaluation. \url{https://www.amazon.com/gp/product/1946728039/}}
\keyword{ package }

