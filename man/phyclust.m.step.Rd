\name{phyclust.m.step}
\alias{phyclust.m.step}
\title{One M-step of phyclust}
\description{
  This is a single M-step of \code{phyclust}, usually following or followed
  by the other E-step.
}
\usage{
phyclust.m.step(X, ret.phyclust = NULL, K = NULL,
    Z.normalized = NULL, substitution.model = NULL,
    identifier = NULL, code.type = NULL)
}
\arguments{
  \item{X}{nid/sid matrix with \eqn{N}{N} rows/sequences and \eqn{L}{L} columns/sites.}
  \item{ret.phyclust}{an object with the class \code{phyclust}.}
  \item{K}{number of clusters.}
  \item{Z.normalized}{posterior probabilities obtained from an E-step.}
  \item{substitution.model}{substitution model.}
  \item{identifier}{identifier.}
  \item{code.type}{code type.}
}
\details{
  \code{X} should be a numerical matrix containing sequence data that
  can be transfered by \code{code2nid} or \code{code2sid}.

  \strong{Note}: gaps \code{-} are not supported yet, drop them from data.

  Either input \code{ret.phyclust} or all other arguments for this function.
  \code{ret.phyclust} can be obtain either from an EM iteration of
  \code{phyclust} or from a E step of \code{phyclust.e.step}.
}
\value{
  This function return an object with class \code{phyclust}.
}
\section{ToDo(s)}{
  \itemize{
    \item implement codes for gaps \code{-}.
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
  \code{\link{phyclust}},
  \code{\link{phyclust.em.step}}.
  \code{\link{phyclust.e.step}}.
}
\examples{
\dontrun{
set.seed(1234)
EMC.1 <- .EMC
EMC.1$EM.iter <- 1
# the same as EMC.1 <- .EMControl(EM.iter = 1)

ret.1 <- phyclust(seq.data.toy$org, 2, EMC = EMC.1)
ret.2 <- phyclust.e.step(seq.data.toy$org, ret.phyclust = ret.1)
ret.3 <- phyclust.m.step(seq.data.toy$org, ret.phyclust = ret.2)
str(ret.3)
}
}
\keyword{programming}
