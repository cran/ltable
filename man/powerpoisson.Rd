\name{powerpoisson}
\alias{PowerPoisson}
\title{Function \emph{PowerPoisson}}
\description{
Performs log-linear and power analyses for constructed tabulated data based on GSL IWLS algorithms
}
\usage{
PowerPoisson(formula, data, scale_min=1, scale_max=5,
effect, p_alpha=0.05, contrasts=NULL )
}
\arguments{
\item{formula}{a symbolic description of the model to be fit.}
\item{data}{name of the data set; object of \emph{data.frame} class}
\item{scale_min}{the smallest number of sample size scale range, 1 signifies the given data sample size (observed total counts).}
\item{scale_max}{the max number of sample size considered in power analysis. 5 by default means 5 times augmented observed counts}
\item{effect}{quoted effect tested by hypothesis; it should be one from the model formula, of second or higher order, introduced by * delimiter, i.e., "y*x", "y1*y2*x1*x2", etc.}
\item{p_alpha}{serves to signify Z to check simulated z-scores against in power analysis, 0.05 by default}
\item{contrasts}{serves to choose types of contrasts to study effects of factors, same with glm), orthogonal polynomials by default}

}

\details{
\itemize{
\item{Performs log-linear modelling with supplied data by using GSL IWLS algorithms.}
\item{Performs power analysis in a given range of sample sizes (scale_min - scale_max).}
\item{Creates object of S4 class \emph{PowerClass}  with accessing methods}
	}
	}
\value{returns object of S4 class \emph{PowerClass} }

\author{
{Ocheredko Oleksandr}
\email{Ocheredko@yahoo.com}
}

\note{
The inspected sample size range defined by scale_min - scale_max automatically is divided into 11 consecutive values investigated by function. Given the results one can change sample size range, for example to scrutinize some particular interval to ensure power and p-value.\cr  Function provides better conditioned variance matrix estimates against function \emph{stats::glm} by the auspices of IWLS GSL algorithms, which is particular important for high order effects and power analysis. Particularly suggestive is to check the model fit first. Jacobian reciprocal condition number near zero indicates solution instability. If \emph{chisq/dof} >> 1,  the error estimates obtained from the covariance matrix will be too small and should be multiplied by square root of \emph{chisq/dof}. Poor fit will result from the use of an inappropriate model and jeopardizes the validity of power analysis. 
}


%% ~Make other sections like Warning with \section{Warning }{....} ~


\seealso{
\code{\link[stats:glm]{glm}}
\code{\link[MCMCglmm:MCMCglmm]{MCMCglmm}}


}

\examples{
require(ltable)
data(tdata, package="ltable")
pres<-PowerPoisson(Counts~smoker +contraceptive +tromb +
contraceptive*tromb, scale_max=1.5, effect="contraceptive*tromb", data=tdata)
print(pres, "model")
print(pres)
plot(pres, stencil=3)
plot(pres)

}


\keyword{utilities}

