\name{table_f}
\alias{table_f}
\title{Function \emph{table_f}}
\description{
Constructs tables of counts and proportions out of data sets.
}
\usage{
table_f(data, datavars, type=1, digits=2, extended=FALSE, MV=FALSE, cb=FALSE)}
\arguments{
\item{data}{name of the data set; object of \emph{data.frame} class}
\item{datavars}{the character string that lists field names separated by comma in the order of presentation in the table: first has its sorted levels rolled out vertically leftmost, the last has its sorted levels spred by columns}
\item{type}{the type of table: 1 (default) - count table; 2 - proportions by rows; 3 - proportions by columns; 4 - frequencies}
\item{digits}{formats output digits number, applied only to proportions, default is 2}
\item{extended}{"+" adds margins of counts, applied only for proportions and frequencies, "-" by default}
\item{MV}{includes missing values into tabulation, operates with type=1 only, FALSE by default}
\item{cb}{TRUE permits to copy the table to clipboard,  FALSE by default}

}

\details{
\itemize{
\item{You can construct table with data set fields of factor, character, logical, and numeric classes.}
\item{To insert table into Word document first open Excell, choose left high corner of placement by mouse click and use Ctrl+V combination or click on the Paste icon (the clipboard), then use Ctrl+C, open Word document, use Ctrl+V to place the table.}
\item{If You want to use clipboard to insert table into Word document use option cb=TRUE.}

	}
	}
\value{returns object of class \emph{data.frame} }

\author{
{Ocheredko Oleksandr}
\email{Ocheredko@yahoo.com}
}

\note{
Abstain from putting continuous variables or too many factor variables into \emph{datavars} list to keep table legible. Put factor variable with numerous levels at the end of the list. 
}


%% ~Make other sections like Warning with \section{Warning }{....} ~


\seealso{
\code{\link[etable:etable]{etable}}
\code{\link[xtable:xtable]{xtable}}
\code{\link[tables:tabular]{tabular}}

}

\examples{
data(sdata, package="ltable")
table_f(sdata, "a")
table_f(sdata, "a", MV=TRUE, extended=TRUE)
table_f(sdata, "a,b,c")
table_f(sdata, "a,b,c,d", type=2, digits=3)
table_f(sdata, "b,c,a,d", MV=TRUE, extended=TRUE, cb=TRUE)
}


\keyword{utilities}

