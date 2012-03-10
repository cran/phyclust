### This file contains functions to call phyclust_se_update in C.

### EM step.
phyclust.se.update <- function(X, se.model = .EMC$se.model,
    se.constant = .EMC$se.constant, ret.phyclust = NULL,
    K = NULL, Eta = NULL, Mu = NULL, pi = NULL, kappa = NULL, Tt = NULL,
    substitution.model = NULL, identifier = NULL, code.type = NULL,
    label = NULL, byrow = TRUE){
  if(is.null(ret.phyclust)){
    if(is.null(K) || is.null(Eta) || is.null(Mu) || is.null(Tt) ||
       is.null(substitution.model) ||
       is.null(identifier) || is.null(code.type)){
      stop("The parameters are not specified correctly.")
    } else{
      ret.phyclust <- list(K = K, Eta = Eta, Mu = Mu,
                           QA = list(pi = pi, kappa = kappa, Tt = Tt,
                                     identifier = identifier),
                           substitution.model = substitution.model,
                           code.type = code.type)
    }
  } else{
    if(class(ret.phyclust) != "phyclust"){
      stop("The ret.phyclust should be in a phyclust class.")
    }
  }

  if(ret.phyclust$code.type != "NUCLEOTIDE"){
    stop("The sequencing error model only supports NUCLEOTIDE data.")
  }

  if(byrow){
    X <- t(X)
  }
  N.X.org <- ncol(X)
  L <- nrow(X)
  K <- ret.phyclust$K

  vect <- convert.QA.to.vect(ret.phyclust)
  label <- check.label(label, N.X.org, K, byrow)

  ret <- .Call("R_phyclust_se_update",
               as.integer(N.X.org),
               as.integer(L),
               as.integer(X),
               as.integer(K),
               as.double(ret.phyclust$Eta),
               as.integer(t(ret.phyclust$Mu)),
               as.double(vect),
               as.integer(which(ret.phyclust$substitution.model ==
                                as.character(.substitution$model)) - 1),
               as.integer(which(ret.phyclust$QA$identifier == .identifier) - 1),
               as.integer(which(ret.phyclust$code.type == .code.type) - 1),
               label,
               as.integer(which(se.model == .se.model) - 1),
               as.double(se.constant),
               PACKAGE = "phyclust")

  if(!is.finite(ret$logL)){
    stop("The logL is not finite.\n")
  }

#  ret$Z.normalized <- ret$bic <- ret$aic <- ret$icl <-
#    ret$class.id <- ret$n.class <- NULL
  ret$class.id <- ret$class.id + 1
  ret$substitution.model <- ret.phyclust$substitution.model
  ret$QA$identifier <- ret.phyclust$QA$identifier
  ret$code.type <- ret.phyclust$code.type

  ret <- translate.ret.se(ret)
  class(ret) <- "phyclust"
  ret
} # End of phyclust.se.update().

