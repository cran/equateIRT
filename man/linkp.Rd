\name{linkp}
\alias{linkp}

\title{
	Linkage Plan
}
\description{
	Calculates the number of common items between a list of forms.
}
\usage{
linkp(coef)
}
\arguments{
	\item{coef}{list of matrices (one for each form) containing item
		parameter estimates. The names of the rows of each matrix should 
		be the names of the items.
	}
}
\value{
	A matrix whose elements indicate the number of common items between
	the forms. On the diagonal of the matrix there are the number of items
	of each form.
}
\author{
	Michela Battauz
}
\seealso{
	\code{\link{est2pl}}, \code{\link{est3pl}}, \code{\link{estrasch}}
}
\examples{
data(est3pl)
linkp(coef = est3pl$coef)
}

