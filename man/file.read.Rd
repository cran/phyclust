\name{file.read}
\alias{read.fasta}
\alias{read.fasta.nucleotide}
\alias{read.phylip}
\alias{read.phylip.nucleotide}
\alias{read.phylip.snp}
\title{Read data from files by formats and return a seq.data}
\description{
  Two major file formats are supported in \pkg{phyclust},
  \code{Format phylip} and \code{Format fasta}.
  \strong{These functions only read files in basic syntax}, and
  return an object in \code{Class seq.data}.
}
\usage{
read.fasta(filename, byrow = TRUE, code.type = .code.type[1])
read.fasta.nucleotide(filename, byrow = TRUE)

read.phylip(filename, byrow = TRUE, code.type = .code.type[1])
read.phylip.nucleotide(filename, byrow = TRUE)
read.phylip.snp(filename, byrow = TRUE)
}
\arguments{
  \item{filename}{a file name where data is read from.}
  \item{byrow}{advanced option, default = TRUE.}
  \item{code.type}{either "NUCLEOTIDE" (default) or "SNP".}
}
\details{
  \code{read.fasta} and \code{read.phylip} are for nucleotide data, and
  \code{read.phylip.snp} is for SNP data.

  \code{byrow} indicates the data will be store by row or not. Usually,
  the default is TRUE. The FALSE is only for advance users
  with careful manipulations and for speeding up the bootstraps.
}
\value{
  Return an object in \code{Class seq.data}.
}
\section{ToDo(s)}{
  \itemize{
    \item implement \code{read.fasta.snp}.
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
  \code{\link{write.fasta}},
  \code{\link{write.phylip}}.
}
\examples{
\dontrun{
# PHYLIP
data.path <- paste(.libPaths()[1], "/phyclust/data/crohn.phy", sep = "")
(my.snp <- read.phylip(data.path, code.type = "SNP"))

# FASTA
data.path <- paste(.libPaths()[1], "/phyclust/data/pony625.fas", sep = "")
(my.pony <- read.fasta(data.path))
}
}
\keyword{programming}
