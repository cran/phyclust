\name{phyclust.Pt}
\alias{phyclust.Pt}
\title{Transition probabilities of phyclust given time}
\description{
  This computes transition probabilities of \code{phyclust} given time.
}
\usage{
phyclust.Pt(Q, Tt, substitution.model = .substitution$model[1],
            log = FALSE)
}
\arguments{
  \item{Q}{a list according to the substitution model.}
  \item{Tt}{total evolution time, \eqn{t}{t}.}
  \item{substitution.model}{substitution model.}
  \item{log}{TRUE for log of transition probabilities.}
}
\details{
  The major models for \code{Q} are:
    \tabular{lll}{
      Model \tab Author and Publication            \tab Parameter \cr
      JC69  \tab Jukes and Cantor 1969.            \tab \eqn{t}{t} \cr
      K80   \tab Kimura 1980.                      \tab \eqn{\kappa, t}{kappa, t} \cr
      F81   \tab Felsenstein 1981.                 \tab \eqn{\pi, t}{pi, t} \cr 
      HKY85 \tab Hasegawa, Kishino, and Yano 1985. \tab \eqn{\pi, \kappa, t}{kappa, pi, t} \cr
    }

  A list of \code{Q} should contains \code{pi}, \code{kappa} based on
  substitution models and code types. \code{Tt} may be separately stored.
  Depending on identifiers, \code{Q}s can be composite to a \code{QA},
  Q matrix array.
}
\value{
  This function returns a transition probability matrix of \code{phyclust}.
}
\section{ToDo(s)}{
  \itemize{
    \item vectorize \code{Tt} for repeated computation in C.
    \item make a general class for \code{Q} and \code{QA}.
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
  \code{\link{.substitution}},
  \code{\link{phyclust}},
  \code{\link{phyclust.em.step}}.
}
\examples{
\dontrun{
Tt <- 0.5

Q <- list(pi = c(0.25, 0.25, 0.25, 0.25), kappa = 0.5)
phyclust.Pt(Q, Tt, "HKY85")

Q <- list(pi = c(0.5, 0.5), kappa = 0.5)
phyclust.Pt(Q, Tt, "SNP_JC69")
}
}
\keyword{programming}
