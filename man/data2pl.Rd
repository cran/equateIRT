\name{data2pl}
\alias{data2pl}
\docType{data}
\title{
	Simulated Data Sets
}
\description{
	Five simulated data sets from a two-parameter logistic model.
}
\usage{
data(data2pl)
}
\format{
	A list of length 5, containing 5 data frames with 5000 dichotomous responses to 20 items.
}
\author{
	Michela Battauz
}

\seealso{
	\code{\link{est2pl}}, \code{\link{import.ltm}}
}
\examples{
data(data2pl)
data2pl[[1]][1:3,]
}
\keyword{datasets}
