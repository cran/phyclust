\name{phyclust.em.step}
\alias{phyclust.em.step}
\title{One EM-step of phyclust}
\description{
  This is a single EM-step of \code{phyclust}.
}
\usage{
phyclust.em.step(X, ret.phyclust = NULL, K = NULL, Eta = NULL,
    Mu = NULL, pi = NULL, kappa = NULL, Tt = NULL,
    substitution.model = NULL, identifier = NULL, code.type = NULL)
}
\arguments{
  \item{X}{nid/sid matrix with \eqn{N}{N} rows/sequences and \eqn{L}{L} columns/sites.}
  \item{ret.phyclust}{an object with the class \code{phyclust}.}
  \item{K}{number of clusters.}
  \item{Eta}{proportion of subpopulations, \eqn{\eta_k}{eta_k}, length = \code{K}, sum to 1.}
  \item{Mu}{centers of subpopulations, dim = \eqn{K\times L}{K*L}, each row is a center.}
  \item{pi}{equilibrium probabilities, each row sums to 1.}
  \item{kappa}{transition and transversion bias.}
  \item{Tt}{total evolution time, \eqn{t}{t}.}
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
  \code{phyclust} or from a M step of \code{phyclust.m.step}.
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
  \code{\link{phyclust.m.step}}.
}
\examples{
\dontrun{
set.seed(1234)
EMC.1 <- .EMC
EMC.1$EM.iter <- 1
# the same as EMC.1 <- .EMControl(EM.iter = 1)

ret.1 <- phyclust(seq.data.toy$org, 2, EMC = EMC.1)
ret.2 <- phyclust.em.step(seq.data.toy$org, ret.phyclust = ret.1)
str(ret.2)
}
}
\keyword{programming}
