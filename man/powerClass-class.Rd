\name{powerClass-class}
\Rdversion{1.1}
\docType{class}
\alias{powerClass-class}

\title{Class \code{"powerClass"}}
\description{
Objects of S4 class \emph{powerClass} are exceptionally suitable for suggested approach to power analysis. Class serves a purpose of container of odds and ends of magnitude of information both on log-linear estimates and fit statistics as well as on the power analysis results, i.e., alpha and beta errors distributions across 11 sample sizes. Class also supported by getters and setters, text and graphic outputs. 
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("powerClass", ...)}.
}
\section{Slots}{
  \describe{
    \item{\code{varnames}:}{Vector of mode \code{"character"} lists names of columns in design matrix.}
    \item{\code{effectsname}:}{Vector of mode \code{"character"} lists names of columns in design matrix that constitute effect under study. Latter is given by arg \emph{effect} in function \emph{PowerPoisson.}}
    \item{\code{cal}:}{Object of class \code{"call"} saves the function call.}
    \item{\code{Ntotal}:}{Vector of mode \code{"numeric"}. Contains sample size of the data, scale_min, scale_max values}
    \item{\code{estim}:}{Object of class \code{"list"} List of 11 lists of log-linear parameters estimates and model fit statistics across 11 sample sizes}
    \item{\code{power1}:}{Object of class \code{"list"}. Contains lists for each column (contrast) of design matrix involved in effect under study. Each such list containes numeric vectors of values of simulated reg.coefficients, z-scores, power. Slot \code{power1} keeps the data pertaining to smallest sample size}
    \item{\code{power2}:}{power2:power11 slots envelop the same structured information across consecutive sample sizes 2:11(largest).}
    \item{\code{power3}:}{-//-}
    \item{\code{power4}:}{-//-}
    \item{\code{power5}:}{-//-}
    \item{\code{power6}:}{-//-}
    \item{\code{power7}:}{-//-}
    \item{\code{power8}:}{-//-}
    \item{\code{power9}:}{-//-}
    \item{\code{power10}:}{-//-}
    \item{\code{power11}:}{-//-}   
    
    }}
\section{Methods}{
  \describe{
    \item{[}{\code{signature(x = "powerClass", i = "character", j = "integer", drop = "logical")}: getter, see Method for Function \code{[}}
    \item{[<-}{\code{signature(x = "powerClass", i = "character", j = "integer", value)}: setter, see Method for Function \code{[<-}}
    \item{plot}{\code{signature(x = "powerClass")}: plots images of z-score and power distributions along the range of sample sizes}
    \item{print}{\code{signature(x = "powerClass")}:  prints estimated log-linear model parameters and fit statistics as well as results of power analysis along the range of sample sizes}
	 }
}
\references{
Ocheredko O.M. MCMC Bootstrap Based Approach to Power and Sample Size Evaluation. \url{https://www.amazon.com/gp/product/1946728039/}}

\author{
{Ocheredko Oleksandr}
\email{Ocheredko@yahoo.com}
}

\examples{
require(ltable)
showClass("powerClass")
new("powerClass")
data(tdata, package="ltable")
## For better illustration You should increase draw and burnin pars
pres<-MCPower(Counts~smoker +contraceptive +tromb +  
contraceptive*tromb, scale_max=1.5, effect="contraceptive*tromb", 
draw=1000, burnin=300, data=tdata)
print(pres)
plot(pres,3)
pres["estim", 1]$betas
pres["power11", 1]$power
pres["power1", 1]$z
}
\keyword{classes}
