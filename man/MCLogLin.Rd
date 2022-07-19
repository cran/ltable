\name{MCLogLin}
\alias{MCLogLin}
\title{Function \emph{MCLogLin}}
\description{
Performs log-linear analyses for constructed tabulated data based on Gibbs sampler with NB2 posterior marginal distribution for counts
}
\usage{
MCLogLin(formula, data, contrasts=NULL, XLB=-100, XUB=100, a=0.1, b=0.1,
DIC=FALSE, pcov=FALSE, draw=10000, burnin=3000 )
}
\arguments{
\item{formula}{a symbolic description of the model to be fit.}
\item{data}{name of the data set; object of \emph{data.frame} class}
\item{contrasts}{serves to choose types of contrasts to study effects of factors, the same with \emph{glm\{stats\}}), orthogonal polynomials by default}
\item{XLB}{the vector of smallest possible values of regression effects \emph{betas}; can be number if pertains to all \emph{betas}.}
\item{XUB}{the vector of largest possible values of regression effects \emph{betas}; can be number if pertains to all \emph{betas}.}
\item{a}{the value of shape parameter of gamma distributed inverce dispersion parameter (\emph{phi}), i.e., \emph{phi~Ga(a,b)}, so that mean(phi)=a/b and var(phi)=a/b^2}
\item{b}{the value of rate (1/scale) parameter of gamma distributed inverce dispersion parameter (\emph{phi}), i.e., \emph{phi~Ga(a,b)}, so that mean(phi)=a/b and var(phi)=a/b^2}
\item{DIC}{requests print of deviance information criteria and its components}
\item{pcov}{requests print of covariance and correlation matricies of the model parameters}
\item{draw}{indicates requested number of samples}
\item{burnin}{indicates requested number of initial samples to discard}

}

\details{
\itemize{
\item{Performs log-linear modelling with supplied data by using Gibbs sampler.}
\item{Printing output includes standard table of parameters estimates, goodness of fit indicators, analysis of residuals. On the prompt it prints the deviance information criteria with its components as well as covariance and correlation matricies of the model parameters }
	}
	}
\value{returns a matrix with columns of chains of sampled values of model parameters (expected counts, regression coefficients, inverce dispersion parameter) to be studied by MCMC facilitating packages (e.g., coda, mcmc, mcmcplot, etc.)}

\author{
{Ocheredko Oleksandr}
\email{Ocheredko@yahoo.com}
}

\note{
Function provides better conditioned variance matrix estimates against function \emph{stats::glm}, which is particular important for high order effects and power analysis. Particularly suggestive is to check the model fit first. Jacobian reciprocal condition number near zero indicates solution instability. If \emph{chisq/n} >> 1,  the error estimates obtained from the covariance matrix will be too small and should be multiplied by square root of \emph{chisq/dof}. Poor fit will result from the use of an inappropriate model and jeopardizes the validity of power analysis.\cr It's recommended to keep difference between pars \emph{draw} and \emph{burnin} at least 3000.   
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
res1<-MCLogLin(Counts~smoker +contraceptive +tromb +
contraceptive*tromb, data=tdata, draw=1500, burnin=500)


data(iris)
iriscut<-with(iris, data.frame(PL=cut(Petal.Length,3),
                               PW=cut(Petal.Width,3)))
irist<-table_f(iriscut,"PL,PW")
irisd<-tableToData(irist, ordered="PL,PW")
res2<-ltable::MCLogLin(Counts~PW+PL+PW*PL, DIC=TRUE, data=irisd,
draw=1500, burnin=500)



}


\keyword{utilities}

