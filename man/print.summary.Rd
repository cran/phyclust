\name{print.object}
\alias{print.ms}
\alias{print.phyclust}
\alias{print.Pt}
\alias{print.RRand}
\alias{print.seq.data}
\alias{print.seqgen}
\alias{summary.phyclust}
\title{Functions for printing or summarizing objects according to classes}
\description{
  Several classes are declared in \pkg{phyclust}, and these are functions
  to print and summary objects.
}
\usage{
print.ms(x, ...)
print.phyclust(x, digits = max(4, getOption("digits") - 3), ...)
print.Pt(x, ...)
print.RRand(x, digits = max(4, getOption("digits") - 3), ...)
print.seq.data(x, ...)
print.seqgen(x, ...)

summary.phyclust(object, ...)
}
\arguments{
  \item{x}{an object with the class attributes.}
  \item{digits}{for printing out numbers.}
  \item{object}{an object with the class attributes.}
  \item{\dots}{other possible options.}
}
\details{
  These are useful functions for summarizing and debugging.

  For \code{ms} and \code{seqgen}, it will show the result as
  standalone versions on the stdout out with line by line.

  For other functions, they only show summaries of objects. Use
  \code{names} or \code{str} to explore the details.
}
\value{
  The results will cat or print on the stdout by default.
}
\references{
  Phylogenetic Clustering Website:
  \url{http://thirteen-01.stat.iastate.edu/snoweye/phyclust/}
}
\author{
  Wei-Chen Chen \email{phyclust@gmail.com}
}
\seealso{
  \code{\link{ms}},
  \code{\link{phyclust}},
  \code{\link{phyclust.Pt}},
  \code{\link{RRand}},
  \code{\link{seqgen}}.
}
\examples{
\dontrun{
# Functions applied by directly type the names of objects.
}
}
\keyword{programming}
