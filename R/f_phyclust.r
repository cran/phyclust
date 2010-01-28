# This file contains functions to call phyclust in C.

phyclust <- function(X, K, EMC = .EMC, manual.id = NULL, byrow = TRUE){
  if(K <= 0){
    stop(K > 0)
  }
  if(byrow){
    X <- t(X)
  }

  N.X.org <- ncol(X)
  L <- nrow(X)

  if(! is.null(manual.id)){
    if(max(manual.id) != K){
      stop("manual.id is not correct.")
    } else{
      manual.id <- manual.id - 1
    }
  } else{
    if(EMC$init.method == "manualMu"){
      stop("manual.id is missing.")
    }
  }

  EMC <- check.EMC(EMC)
  EMC <- translate.EMC(EMC)

  ret <- .Call("R_phyclust",
               as.integer(N.X.org),
               as.integer(L),
               as.integer(K),
               as.integer(X),
               EMC,
               as.integer(manual.id),
               PACKAGE = "phyclust")

  if(!is.finite(ret$logL)){
    stop("logL is not finite.\n")
  }

  ret$class.id <- ret$class.id + 1
  ret <- translate.ret(ret, EMC)
  class(ret) <- "phyclust"
  ret
} # End of phyclust().


### For internal used.
check.EMC <- function(EMC){
  if(EMC$init.method %in% c("NJ", "PAM", "manualMu") &&
     EMC$init.procedure != "exhaustEM"){
    my.cat("init procedure: ", EMC$init.procedure, " -> exhaustEM",
           " (", EMC$init.method, ")\n")
    EMC$init.procedure <- "exhaustEM"
    EMC$exhaust.iter <- 1
  }

  if((EMC$code.type == "SNP") &&
     (!EMC$substitution.model %in% c("SNP_JC69", "SNP_F81", "E_SNP_F81"))){
    my.cat("substitution model: ", EMC$substitution.model, " -> SNP_JC69",
           " (", EMC$code.type, ")\n")
    EMC$substitution.model <- "SNP_JC69"
  }

  if(EMC$substitution.model %in% c("SNP_JC69", "SNP_F81", "E_SNP_F81") &&
     EMC$edist.model != "D_HAMMING"){
    my.cat("edist model: ", EMC$edist.model, " -> D_HAMMING",
           " (", EMC$substitution.model, ")\n")
    EMC$edist.model <- "D_HAMMING"
  }

  EMC
} # End of check.EMC().

translate.ret <- function(ret, EMC = NULL){
  ret$Mu <- matrix(ret$Mu, nrow = ret$K, byrow = TRUE)
  if(! is.null(ret$Z.normalized)){
    ret$Z.normalized <- matrix(ret$Z.normalized, nrow = ret$N.X.org, byrow = TRUE)
  }
  ret$QA$pi <- matrix(ret$QA$pi, nrow = ret$K, byrow = TRUE)
  ret$QA$kappa <- matrix(ret$QA$kappa, nrow = ret$K, byrow = TRUE)
  ret$QA$Tt <- matrix(ret$QA$Tt, nrow = ret$K, byrow = TRUE)

  if(!is.null(EMC)){
    ret$init.procedure <- .init.procedure[EMC$init.procedure + 1]
    ret$init.method <- .init.method[EMC$init.method + 1]
    ret$substitution.model <-
      as.character(.substitution$model[EMC$substitution.model + 1])
    ret$edist.model <- .edist.model[EMC$edist.model + 1]
    ret$QA$identifier <- .identifier[EMC$identifier + 1]
    ret$code.type <- .code.type[EMC$code.type + 1]
    ret$em.method <- .em.method[EMC$em.method + 1]
    ret$boundary.method <- .boundary.method[EMC$boundary.method + 1]
  }

  if(ret$QA$identifier == "EE"){
    ret$QA$pi <- matrix(ret$QA$pi[1,], nrow = 1)
    ret$QA$kappa <- ret$QA$kappa[1]
    ret$QA$Tt <- ret$QA$Tt[1]
  } else if(ret$QA$identifier == "EV"){
    ret$QA$pi <- matrix(ret$QA$pi[1,], nrow = 1)
    ret$QA$kappa <- ret$QA$kappa[1]
  } else if(ret$QA$identifier == "VE"){
    ret$QA$Tt <- ret$QA$Tt[1]
  }

  if(ret$code.type == "NUCLEOTIDE"){
    colnames(ret$QA$pi) <- as.character(.nucleotide$code[1:ncol(ret$QA$pi)])
  } else if(ret$code.type == "SNP"){
    colnames(ret$QA$pi) <- as.character(.snp$code[1:ncol(ret$QA$pi)])
  }
  rownames(ret$QA$pi) <- paste("k=", 1:nrow(ret$QA$pi), sep = "")

  if(ret$substitution.model %in% c("JC69", "K80", "SNP_JC69")){
    ret$QA$pi <- NULL
  }
  if(ret$substitution.model %in%
     c("JC69", "F81", "SNP_JC69", "SNP_F81", "E_F81", "E_SNP_F81")){
    ret$QA$kappa <- NULL
  }

  ret
} # End of translate.ret().

translate.EMC <- function(EMC){
  EMC$exhaust.iter <- as.integer(EMC$exhaust.iter)
  EMC$fixed.iter <- as.integer(EMC$fixed.iter)
  EMC$short.iter <- as.integer(EMC$short.iter)
  EMC$EM.iter <- as.integer(EMC$EM.iter)
  EMC$short.eps <- as.double(EMC$short.eps)
  EMC$EM.eps <- as.double(EMC$EM.eps)

  EMC$cm.reltol <- as.double(EMC$cm.reltol)
  EMC$cm.maxit <- as.integer(EMC$cm.maxit)

  EMC$nm.abstol.Mu.given.QA <- as.double(EMC$nm.abstol.Mu.given.QA)
  EMC$nm.abstol.QA.given.Mu <- as.double(EMC$nm.abstol.QA.given.Mu)
  EMC$nm.reltol.Mu.given.QA <- as.double(EMC$nm.reltol.Mu.given.QA)
  EMC$nm.reltol.QA.given.Mu <- as.double(EMC$nm.reltol.QA.given.Mu)
  EMC$nm.maxit.Mu.given.QA <- as.integer(EMC$nm.maxit.Mu.given.QA)
  EMC$nm.maxit.QA.given.Mu <- as.integer(EMC$nm.maxit.QA.given.Mu)
  EMC$est.non.seg.site <- as.integer(EMC$est.non.seg.site)

  EMC$min.n.class <- as.integer(EMC$min.n.class)

  if(EMC$init.procedure[1] %in% .init.procedure){
    EMC$init.procedure <- as.integer(which(EMC$init.procedure ==
                                           .init.procedure) - 1)
  } else{
    stop("The initial procedure is not found.")
  }

  if(EMC$init.method[1] %in% .init.method){
    EMC$init.method <- as.integer(which(EMC$init.method == .init.method) - 1)
  } else{
    stop("The initial method is not found.")
  }

  if(EMC$substitution.model[1] %in% as.character(.substitution$model)){
    EMC$substitution.model <-
      as.integer(which(EMC$substitution.model ==
                       as.character(.substitution$model)) - 1)
  } else{
    stop("The substitution model is not found.")
  }

  if(EMC$edist.model[1] %in% .edist.model){
    EMC$edist.model <- as.integer(which(EMC$edist.model == .edist.model) - 1)
  } else{
    stop("The edist model is not found.")
  }

  if(EMC$identifier[1] %in% .identifier){
    EMC$identifier <- as.integer(which(EMC$identifier == .identifier) - 1)
  } else{
    stop("The identifier is not found.")
  }

  if(EMC$code.type[1] %in% .code.type){
    EMC$code.type <- as.integer(which(EMC$code.type == .code.type) - 1)
  } else{
    stop("The code type is not found.")
  }

  if(EMC$em.method[1] %in% .em.method){
    EMC$em.method <- as.integer(which(EMC$em.method == .em.method) - 1)
  } else{
    stop("The em method is not found.")
  }

  if(EMC$boundary.method[1] %in% .boundary.method){
    EMC$boundary.method <- as.integer(which(EMC$boundary.method ==
                                            .boundary.method) - 1)
  } else{
    stop("The boundary method is not found.")
  }

  EMC$max.init.iter <- as.integer(EMC$max.init.iter)

  EMC
} # End of translate.EMC().

