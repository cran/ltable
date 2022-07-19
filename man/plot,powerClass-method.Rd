\name{plot,powerClass-method}
\docType{methods}
\alias{plot,powerClass-method}
\title{ Method for Function \code{plot}}
\description{
 Method for function \code{plot} with \cr
 \code{signature(x = "powerClass")}
}
\section{Methods}{
\describe{
\item{\code{signature(x = "powerClass")}}{
Method for function \code{plot} for object of S4 class \emph{powerClass}.\cr

}
}}
\usage{
\S4method{plot}{powerClass}(x, stencil, \dots)
}


\arguments{
\item{x}{the name of \emph{powerClass} object.}
\item{stencil}{an optional arg containing 4 choices of print: missing(default), 1, 2, 3. See details.}
\item{\dots}{not used}
}
\details{
The second argument \emph{stencil} controls "what and how" to plot. \emph{stencil=missing} (default) plots stand-alone images of z-score and power distributions along the range of sample sizes (see \emph{print-method} for details on the range).\cr
\emph{stencil=1} chooses z-score distributions to plot in stand-along fashion.\cr
\emph{stencil=2} chooses power distributions to plot in stand-along fashion.\cr
\emph{stencil=3} controls to plot z-score and power distributions paired alongside.\cr
Also, Q0.05, Q01, Q0.5 (median) quantiles are graphed in lines.   
}

\examples{
require(ltable)
data(tdata, package="ltable")
## For better illustration You should increase draw and burnin pars
pres<-MCPower(Counts~smoker +contraceptive +tromb +  
contraceptive*tromb, scale_max=1.5, effect="contraceptive*tromb", 
data=tdata, draw=1000, burnin=300)
plot(pres)
plot(pres, stencil=3)
}
\keyword{methods}
