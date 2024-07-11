\name{renewal}
\alias{renewal}
\title{Function \emph{renewal}}
\description{
Estimates MCMC chain convergence based on renewal theory.
}
\usage{
renewal(x, StatesNum=10, Astate=NULL, nForStart = 3000, epsilon1=0.05, epsilon2=0.05)
}
\arguments{
\item{x}{name of the numeric vector of the chain; object of \emph{numeric} class}
\item{StatesNum}{positive integer that indicates the number of states to classify chain values into, max value is 100 s.t. values of chain classified into 100 ordered states; to choose StatesNum consider the adequate number of unique intervals chain values can be suitably represented by}
\item{Astate}{positive integer from 1 to StateNum that defines state to base calculation on; there may be small differences in results under different base states, see help for details; default is median state value}
\item{nForStart}{positive integer that indicates the maximum distance from beginning of chain to consider in finding Start value, 3000 by default}
\item{epsilon1}{upper bound of closeness to stationary distribution \eqn{\pi_X}; it is used to calculate the starting position of the chain \emph{nstat}
\deqn{\|P^{nstat}_{x0} - \pi_X \| \leqslant \epsilon_1}}
\item{epsilon2}{ is upper bound of variance of estimator
\deqn{Var(\frac{1}{\ell}\sum^{nvar}_{nstat+1}X_k - E_{\pi_X}) \leqslant \epsilon_2}}
}

\details{
\itemize{
\item{You can ascertaine if the sampled transition probability is close to the determined stationary probability of Markov Chain and how many iterations should be used in order to minimize the error of estimator.}
\item{Algorithm is based on renewal theory and implements the concept of the so called "secondary chain". It is supported by strong mathematical reasoning and the obtained solutions are strict, i.e. they are not asymptotic. Hence, this method is not biased by additional error provided by limit theorems.}
\item{If the calculated chain length surpasses length of supplied chain whatever chain length is provided it is the indication that chain based estimator deviates from stationary distribution based more than indicated by \emph{epsilon2}. }   
\item{The larger the \emph{epsilon1} the larger the starting value of the chain. The larger the \emph{epsilon2} the lengthier is the chain. }
	}
	}
\value{returns two positive integers, starting value of the chain and length of the chain to be used from starting value. }

\author{
{Ocheredko Oleksandr}
\email{Ocheredko@yahoo.com}
}

\note{
Abstain from repeated supply of lengthier chain. It does not correct the problem of nonstationarity.  
}

\examples{
require(ltable)
data(tdata, package="ltable")
res1<-MCLogLin(Counts~smoker +contraceptive +tromb +
contraceptive*tromb, data=tdata, draw=5000, burnin=500)
renewal(res1[,1], Astate=5)
renewal(res1[,14], Astate=1)

}


\keyword{utilities}

