\name{[,powerClass-method}
\docType{methods}
\alias{[,powerClass-method}
\title{ Method for Function \code{[}}
\description{
 Method for function \code{[} with \cr
 \code{signature(x = "powerClass")}
}
\section{Methods}{
\describe{
\item{\code{signature(x = "powerClass", i = "character", j = "integer", drop = "logical")}}{
Method for function \code{[} for object of S4 class \emph{powerClass}.\cr 
}
}
}

\usage{
\S4method{[}{powerClass}(x, i, j, drop)
}


\arguments{
\item{x}{the name of \emph{powerClass} object.}
\item{i}{the name of the slot of the object}
\item{j}{picks up j-th element of the list in slot with name &i.}
\item{drop}{not used}
}

\details{
Method provides access to slots of \emph{powerClass} object. Its structure delivered in \emph{powerClass-class} index. Access to particular vectors of lists supplied with $ operator. For example, log-linear reg.coefficients estimates of smallest size data accessible by obj["estim", 1]$betas, errors can be obtained by analogue: obj["estim", 1]$errors. Power values extraction slightly differs: obj["power11", 1]&power extracts power values vector for 1st effect given 11th (largest) sample size. By analogue we get vector of z-scores for second effect given smallest sample size by obj["power1", 2]&z.      
See \emph{powerClass-class} index.}

\examples{
require(ltable)
require(ltable)
data(tdata, package="ltable")
## For better illustration You should increase draw and burnin pars
pres<-MCPower(Counts~smoker +contraceptive +tromb +  
contraceptive*tromb, scale_max=1.5, effect="contraceptive*tromb", 
data=tdata, draw=1000, burnin=300)
# get call
pres["cal"]
# get effect contrasts names
pres["effectsname"]
# get Jacobian reciprocal condition number in smallest sample
pres["estim",1]$Jacobian_rcnumber
# get chisq/n in smallest sample
pres["estim",1]$chi_sq
# get LogLikelihood
pres["estim",1]$LL
# get initial deviation between observed and expected counts
pres["estim",1]$dev0
# get final deviation between observed and expected counts
pres["estim",1]$dev

}
\keyword{methods}
