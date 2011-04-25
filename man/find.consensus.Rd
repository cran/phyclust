\name{find.consensus}
\alias{find.consensus}
\title{Find the consensus sequence}
\description{
  Based on the input data, this function will seach all data along
  all sites to find a ceonsensus sequence which may be or may not be
  one of the data.
}
\usage{
find.consensus(X, code.type = .code.type[1])
}
\arguments{
  \item{X}{nid/sid matrix with \eqn{N}{N} rows/sequences and \eqn{L}{L} columns/sites.}
  \item{code.type}{either "NUCLEOTIDE" (default) or "SNP".}
}
\details{
  \code{X} should be a numerical matrix containing sequence data that
  can be transfered by \code{code2nid} or \code{code2sid}.
}
\value{
  A vector containing the consensus sequence with length \eqn{L}{L}
  will be returned.
}
\references{
  Phylogenetic Clustering Website:
  \url{http://thirteen-01.stat.iastate.edu/snoweye/phyclust/}
}
\author{
  Wei-Chen Chen \email{phyclust@gmail.com}
}
\seealso{
  \code{\link{plotdots}}.
}
\examples{
\dontrun{
find.consensus(seq.data.toy$org)
}
}
\keyword{programming}
