\name{ltable-package}
\alias{ltable}
\docType{package}
\title{
ltable}

\description{
\bold{Constructs tables of counts and proportions out of data sets. Performs log-linear and power analyzes of tabulated data}
}


\details{

Gibbs sampling based log-linear analysis features some advantages against \emph{glm \{stats\}}, first of all due to fixing overdispersion by NB2 posterior marginal distribution of counts that insures distinctly less biased covariance estimates, pivot issue for implemented power analysis. In some instances hypothesis testing of higher order effects disagrees with that of \emph{glm \{stats\}} on account of larger NB2 model based errors estimates. Another though related enhancement is distinct better fit assessed by sum of squared differences between observed and expected counts. Results of power analysis backed up with MCMC BUGS delivered approach (reference 2).  

}


\author{
{Ocheredko Oleksandr}
\email{Ocheredko@yahoo.com}
}

\note{
You can:
\enumerate{
\item{construct tables with data set fields of factor, character, logical, and numeric classes;}
\item{insert tables into Excel and Word documents using clipboard, into LaTeX, HTML, Markdown and reStructuredText documents by the knitr::kable agency;}
\item{perform Gibbs sampling based log-linear analysis;}
\item{perform power analysis of selected effect.}
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
Ocheredko O.M. MCMC Bootstrap Based Approach to Power and Sample Size Evaluation. \url{https://www.amazon.com/gp/product/1946728039/}}
\keyword{ package }

