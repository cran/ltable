\name{print,powerClass-method}
\docType{methods}
\alias{print,powerClass-method}
\title{ Method for Function \code{print}  }
\description{
 Method for function \code{print} with \cr
 \code{signature(x = "powerClass")}
}
\section{Methods}{
\describe{
\item{\code{signature(x = "powerClass")}}{
Method for function \code{print} for object of S4 class \emph{powerClass}.\cr The second argument \emph{choice} controls information to print. It's advisable to start printing with arg \emph{choice="model"}. Besides estimated log-linear model parameters, fit statistics printed for input data given arg \emph{scale_min}=1 in function \emph{PowerPoisson}. Otherwise, it prints results for augmented scale_min*data counts. Of particular importance is \emph{Jacobian reciprocal condition number} and \emph{chisq/dof}. See details.\cr
Arg \emph{choice="power"} prints results of power analysis in given range of sample size regulated by args \emph{scale_min, scale_max} in function \emph{PowerPoisson}. These are multiplyers for observed data counts. Range is divided into 11 even-spaced subsequent sample sizes. Each is described in printed quantiles (Q0.025, Q0.05, Q0.1, Q0.2, Q0.3, Q0.4, Q0.5) of power and z-score distributions. It's suggestive to use Q0.025 in making decision.  Given the results one can change sample size range, for example to scrutinize some particular interval to ensure power and p-value.    

}
}}
\usage{
\S4method{print}{powerClass}(x, choice, \dots)
}


\arguments{
\item{x}{the name of \emph{powerClass} object.}
\item{choice}{an optional arg containing two choices of print: "power" (by default) prints the results of power analysis, while "model" prints estimated log-linear model parameters and fit statistics.}
\item{\dots}{not used}
}
\details{
Fit statistic \emph{Jacobian reciprocal condition number} measures the inverse sensitivity of the solution to small perturbations in the input data. It tends to zero as \emph{J} tends to singularity indicating solution instability. \cr
The value of ch-squared per degree of freedom \emph{chisq/dof} approximately 1 indicates a good fit. If \emph{chisq/dof} >> 1  the error estimates obtained from the covariance matrix will be too small and should be multiplied by square root of \emph{chisq/dof}. \cr
Poor fit will result from the use of an inappropriate model.\cr 
BEWARE: Poor fit jeopardizes the validity of power analysis.
}
\examples{
require(ltable)
data(tdata, package="ltable")
## For better illustration You should increase draw and burnin pars
pres<-MCPower(Counts~smoker +contraceptive +tromb +
contraceptive*tromb, scale_min = 0.5, scale_max=1.5, 
effect="contraceptive*tromb", data=tdata, draw=1000, burnin=300)
print(pres, "model")
print(pres, "power")

}
\keyword{methods}

