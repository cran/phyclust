\name{phyclust.se.update}
\alias{phyclust.se.update}
\title{Update phyclust Results by the Sequencing Error Model}
\description{
  Since \code{phyclust.se} is difficult to optimize on a constrained high
  dimension parameter space, the \code{phyclust} is relatively easier to
  find a better result, as well as the \code{find.best} function.

  This function will use the \code{phyclust} result as initial parameters and
  perform a sequencing error model. All parameters (Eta, Mu, Q, ...) in this
  function will be updated through the EM algorithm as \code{phyclust.se}.

  Typically, this function run on the \code{find.best} results will yield
  a better result than on the \code{phyclust}.
}
\usage{
phyclust.se.update(X, se.model = .EMC$se.model, se.constant = .EMC$se.constant,
    ret.phyclust = NULL, K = NULL, Eta = NULL,
    Mu = NULL, pi = NULL, kappa = NULL, Tt = NULL,
    substitution.model = NULL, identifier = NULL, code.type = NULL,
    label = NULL, byrow = TRUE)
}
\arguments{
  \item{X}{nid/sid matrix with \eqn{N}{N} rows/sequences and \eqn{L}{L} columns/sites.}
  \item{se.model}{sequencing error model, default = "CONVOLUTION".}
  \item{se.constant}{constrained constant, default = 1e-2.}
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
  \item{label}{label of sequences for semi-supervised clustering.}
  \item{byrow}{advanced option for \code{X}, default = TRUE.}
}
\details{
  All the input arguments are the same as the inputs of the function
  \code{phyclust.em.step}.
}
\value{
  This function return an object with class \code{phyclust}.
}
\references{
  Phylogenetic Clustering Website:
  \url{http://thirteen-01.stat.iastate.edu/snoweye/phyclust/}
}
\author{
  Wei-Chen Chen \email{phyclust@gmail.com}
}
\seealso{
  \code{\link{phyclust.se}},
  \code{\link{phyclust}},
  \code{\link{find.best}}.
}
\examples{
\dontrun{
set.seed(1234)
X <- seq.data.toy$org

(ret.1 <- find.best(X, 4))
(ret.2 <- phyclust.se.update(X, ret.phyclust = ret.1))
(ret.3 <- phyclust.se.update(X, ret.phyclust = ret.1, se.constant = 1e-3))
}
}
\keyword{phyclust}

