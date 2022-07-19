\name{MCPower}
\alias{MCPower}
\title{Function \emph{MCPower}}
\description{
Performs power analyses for constructed tabulated data based on based on Gibbs sampler with NB2 posterior marginal distribution for counts
}
\usage{
MCPower(formula, data, contrasts=NULL, XLB=-100, XUB=100, a=0.1, b=0.1,
scale_min=1, scale_max=5, effect, p_alpha=0.05, draw=10000, burnin=3000)
}
\arguments{
\item{formula}{a symbolic description of the model to be fit.}
\item{data}{name of the data set; object of \emph{data.frame} class}
\item{contrasts}{serves to choose types of contrasts to study effects of factors, same with \emph{glm\{stats\}}), orthogonal polynomials by default}
\item{XLB}{the vector of smallest possible values of regression effects \emph{betas}; can be number if pertains to all \emph{betas}.}
\item{XUB}{the vector of largest possible values of regression effects \emph{betas}; can be number if pertains to all \emph{betas}.}
\item{a}{the value of shape parameter of gamma distributed inverce dispersion parameter (\emph{phi}), i.e., \emph{phi~Ga(a,b)}, so that mean(phi)=a/b and var(phi)=a/b^2}
\item{b}{the value of rate (1/scale) parameter of gamma distributed inverce dispersion parameter (\emph{phi}), i.e., \emph{phi~Ga(a,b)}, so that mean(phi)=a/b and var(phi)=a/b^2}
\item{scale_min}{the smallest number of sample size scale range, 1 signifies the given data sample size (observed total counts).}
\item{scale_max}{the max number of sample size considered in power analysis. 5 by default means 5 times augmented observed counts}
\item{effect}{quoted effect tested by hypothesis; it should be one from the model formula, of second or higher order, introduced by * delimiter, i.e., "y*x", "y1*y2*x1*x2", etc.}
\item{p_alpha}{serves to signify Z to check simulated z-scores against in power analysis, 0.05 by default}
\item{draw}{indicates requested number of samples}
\item{burnin}{indicates requested number of initial samples to discard}

}

\details{
\itemize{
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
The inspected sample size range defined by scale_min - scale_max automatically is divided into 11 consecutive values investigated by function. Given the results one can change sample size range, for example to scrutinize some particular interval to ensure power and p-value.\cr  Function provides better conditioned variance matrix estimates against function \emph{stats::glm} by the auspicity of NB2 dispersion parameter, coping with overdispersion in counts distribution, which is particular important for high order effects and power analysis. Particularly suggestive is to check the model fit first. Jacobian reciprocal condition number near zero indicates solution instability. If \emph{chisq/n} >> 1,  the error estimates obtained from the covariance matrix will be too small and should be multiplied by square root of \emph{chisq/dof}. Poor fit will result from the use of an inappropriate model and jeopardizes the validity of power analysis. \cr The drawback is failure to tackle singularity of order 5 or higher of Hessian matrix. Code returns error "Sorry, can't proceed with singular Hessian matrix." On such rare occasions please use ltable v.2.0.1 available for Unix (MacOS) machines. Function PowerPoisson performs log-linear and power analyses based on Levenberg-Marquardt algorithm which is distribution-free (so, Poisson in name of function is misleading). The only inconvenience is that GSL: GNU Scientific Library has to be installed first. \cr See-saw dynamic of either power or test curves is caused by Jacobian singularity, that indicates solution instability. \cr Flat profiles given low test or power values are indicative for insignificance of tested effect. \cr Flat profiles with z-values above 2 or power values that exceed 0.8 are indicative for significance of tested effect. On such occasions decrease both scale parameters to inspect smaller samples.   
}


%% ~Make other sections like Warning with \section{Warning }{....} ~


\seealso{
\code{\link[stats:glm]{glm}}
\code{\link[MCMCglmm:MCMCglmm]{MCMCglmm}}


}

\examples{
require(ltable)
data(tdata, package="ltable")
## For better illustration You should increase draw and burnin pars
pres1<-MCPower(Counts~smoker +contraceptive +tromb +contraceptive*tromb,
scale_min=0.5, scale_max=1.5, effect="contraceptive*tromb", data=tdata, 
draw=1000, burnin=300)
print(pres1, "model")
print(pres1)
plot(pres1, stencil=3)

}


\keyword{utilities}

