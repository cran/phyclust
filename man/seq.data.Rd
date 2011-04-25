\name{seq.data}
\alias{seq.data.toy}
\alias{seq.data.missing}
\alias{Class seq.data}
\alias{Data Toy}
\docType{data}
\title{A toy dataset in class seq.data}
\description{
   A toy dataset, \code{seq.data.toy}, with 100 nucleotide sequences
   and 200 sites in 4 clusters.
   \code{seq.data.missing} contains some missings indicated by "-".
}
\format{
  This data contains a list with a \code{seq.data} structure described
  in the \sQuote{Details} section.
}
\details{
  A toy dataset is generated to demonstrate \pkg{phyclust}.
  It has 100 nucleotide sequences and 200 sites in 4 clusters
  where the ancestral tree height 0.15 and the descendant tree height 0.09,
  and sequences are evolved by a HKY85 model.

  The structre of class \code{seq.data} is a list containing:
  \tabular{ll}{
    \code{code.type} \tab either "NUCLEOTIDE" or "SNP". \cr
    \code{info}      \tab header for phylip file. \cr
    \code{nseq}      \tab number of sequences, \eqn{N}{N}. \cr
    \code{seqlen}    \tab length of sequences, \eqn{L}{L}. \cr
    \code{seqname}   \tab sequence names. \cr
    \code{org.code}  \tab original codes, dim = \eqn{N\times L}{N*L}. \cr
    \code{org}       \tab transfered ids, dim = \eqn{N\times L}{N*L}. \cr
  }
}
\references{
  Phylogenetic Clustering Website:
  \url{http://thirteen-01.stat.iastate.edu/snoweye/phyclust/}
}
\author{
  Wei-Chen Chen \email{phyclust@gmail.com}
}
\seealso{
  \code{\link{print.seq.data}}.
}
\examples{
\dontrun{
seq.data.toy
seq.data.missing

par(mfrow = c(1, 2))
plotdots(seq.data.toy$org)
plotdots(seq.data.missing$org)
}
}
\keyword{datasets}
