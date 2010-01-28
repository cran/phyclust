\name{standard.code}
\alias{.nucleotide}
\alias{.snp}
\docType{data}
\title{Standard codes and ids for nucleotides and SNPs}
\description{
  Standard codes and ids for nucleotides and SNPs.
  Both objects are used to transfer data.
  \strong{This is a read-only object and the elemental order is followed in C.}
}
\usage{
.nucleotide
.snp
}
\format{
  Both objects are data frames containing ids and codes.
}
\details{
  \strong{Note}: All ids are coding started from \strong{0}.

  Nucleotides, A, G, C, T, and - have codes 0, 1, 2, 3, and 4.
  SNPs, 1, 2, and - have code codes 0, 1, and 2.

  These are two objects are in data frames unlike other internal objects due
  to heavily used in processing data. The original data should be transfered
  to numerical codes in order to be passed to C codes. In C codes, we use
  integers, 0, 1, 2, ... for coding nucleotides or SNPs.

  Data frames use factor formats as the default, and \code{as.character} is
  the way to transfer to the characters.
}
\references{
  Phylogenetic Clustering Website:
  \url{http://thirteen-01.stat.iastate.edu/snoweye/phyclust/}
}
\author{
  Wei-Chen Chen \email{phyclust@gmail.com}
}
\seealso{
  \code{\link{.show.option}},
  \code{\link{code2nid}},
  \code{\link{nid2code}},
  \code{\link{snp2sid}},
  \code{\link{sid2snp}}.
}
\examples{
\dontrun{
.nucleotide
.snp
}
}
\keyword{data}
