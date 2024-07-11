\name{tableToData}
\alias{tableToData}
\title{Function \emph{tableToData}}
\description{
Constructs \emph{data.frames} that fit \emph{glm\{stats\}} or \emph{MCLogLin\{ltable\}}, \emph{MCPower\{ltable\}} modelling out of tables created with function \emph{table_f\{ltable\}}.
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
\item{Check the input and output. It's not a problem to have huge counts together with zero counts for NB2 model used in \emph{ltable}. Still good practice to proceed with power analisis is to have data without zero counts. It's in no way detrimental as in the case of the Poisson GLM, having the mean and variance equality. The implication with Poisson GLM is that as the mean tends to zero, so must the variance. Still we do have some uncertainty about this fitted value. Of the same nature (but converse) problem is with cells of large counts.}
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
data(iris)
iriscut<-with(iris, data.frame(PL=cut(Petal.Length,3),
                               PW=cut(Petal.Width,3)))
irist<-table_f(iriscut,"PL,PW")
irisd<-tableToData(irist, ordered="PL,PW")

}


\keyword{utilities}

