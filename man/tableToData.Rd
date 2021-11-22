\name{tableToData}
\alias{tableToData}
\title{Function \emph{tableToData}}
\description{
Constructs data.frames that fit glm or PowerPoisson modelling out of tables created with function \emph{table_f}.
}
\usage{
tableToData(tname, numerictype="", orderedtype="")
}
\arguments{
\item{tname}{name of the tables created with function \emph{table_f}; object of \emph{data.frame} class}
\item{numerictype}{the character string that lists variable names separated by comma to be transformed to numeric class. Variable "Counts" shouldn't be listed}
\item{orderedtype}{the character string that lists variable names separated by comma to be transformed to ordered factor class. Variable "Counts" shouldn't be listed}
}
\details{
\itemize{
\item{Variables of character and logical classes shape the same design as does the factor class, therefore there is no need to change them to factors.}
\item{Check the input and output. Good practice is to have data without zero counts. In the Poisson GLM, the mean and variance are the same thing. The implication of this is that as the mean tends to zero, so must the variance. Still we do have some uncertainty about this fitted value. Of the same nature (but converse) problem is with cells of large counts.}
\item{You can build the data by hand and skip this functionality.}
	}
	}
\value{returns object of class \emph{data.frame} }

\author{
{Ocheredko Oleksandr}
\email{Ocheredko@yahoo.com}
}


\seealso{
\code{\link[stats:reshape]{reshape}}



}

\examples{
require(ltable)
data(sdata, package="ltable")
stab<-table_f(sdata, "a,b,c")
sdat<-tableToData(stab, orderedtype="c")
res<-PowerPoisson(Counts~a+b+c+a*c, effect="a*c", data=sdat)
plot(res, st=3)
}


\keyword{utilities}

