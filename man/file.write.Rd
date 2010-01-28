\name{file.write}
\alias{write.fasta}
\alias{write.fasta.nucleotide}
\alias{write.phylip}
\alias{write.phylip.nucleotide}
\alias{write.phylip.snp}
\title{Write data to files by formats}
\description{
  Two major file formats are supported in \pkg{phyclust},
  \code{Format phylip} and \code{Format fasta}.
  \strong{These functions only write files in basic syntax.}
}
\usage{
write.fasta(seqdata, filename, classid = NULL, seqname = NULL,
    width.line = 60, lower.case = FALSE, code.type = .code.type[1])
write.fasta.nucleotide(seqdata, filename, classid = NULL,
    seqname = NULL, width.line = 60, lower.case = FALSE)

write.phylip(seqdata, filename, classid = NULL, seqname = NULL,
    width.seqname = 10, width.line = 60, lower.case = FALSE,
    code.type = .code.type[1])
write.phylip.nucleotide(seqdata, filename, classid = NULL,
    seqname = NULL, width.seqname = 10, width.line = 60,
    lower.case = FALSE)
write.phylip.snp(seqdata, filename, classid = NULL,
    seqname = NULL, width.seqname = 10, width.line = 60)
}
\arguments{
  \item{seqdata}{a matrix contains sequence ids as \code{X} in \code{phyclust}.}
  \item{filename}{a file name where data is written to.}
  \item{classid}{class id of sequences.}
  \item{seqname}{sequence names.}
  \item{width.seqname}{number of characters of sequence names to be stored.}
  \item{width.line}{width of lines for breaking lines.}
  \item{lower.case}{use lower case of letters to write}
  \item{code.type}{either "NUCLEOTIDE" (default) or "SNP".}
}
\details{
  \code{write.fasta} and \code{write.phylip} are for nucleotide data, and
  \code{write.phylip.snp} is for SNP data.
}
\value{
  Save a text file.
}
\section{ToDo(s)}{
  \itemize{
    \item implement \code{write.fasta.snp}.
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
  \code{\link{read.fasta}},
  \code{\link{read.phylip}}.
}
\examples{
\dontrun{
# PHYLIP
data.path <- paste(.libPaths()[1], "/phyclust/data/crohn.phy", sep = "")
my.snp <- read.phylip(data.path, code.type = "SNP")
write.phylip(my.snp$org, "new.crohn.phy", code.type = "SNP")

# FASTA
data.path <- paste(.libPaths()[1], "/phyclust/data/pony625.fas", sep = "")
(my.pony <- read.fasta(data.path))
write.phylip(my.pony$org, "new.pony.phy")
}
}
\keyword{programming}
