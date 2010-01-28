### This file contains functions for computing P(t) = e^{Qt}.
# NUCLEOTIDE: Q <- list(pi = c(0.25, 0.25, 0.25, 0.25), kappa = 0.5)
# SNP: Q <- list(pi = c(0.5, 0.5), kappa = 0.5)

phyclust.Pt <- function(Q, Tt, substitution.model = .substitution$model[1],
    log = FALSE){
  sub.model <- which(substitution.model[1] == as.character(.substitution$model))
  code.type <- which(as.character(.substitution$code.type[sub.model]) ==
                                  .code.type)
  ret <- .Call("R_phyclust_logPt",
               as.double(Q$pi),
               as.double(Q$kappa),
               as.double(Tt),
               as.integer(code.type - 1),
               as.integer(sub.model - 1),
               PACKAGE = "phyclust")
  if(!log){
    ret <- exp(ret)
  }
  if(as.character(.substitution$code.type[sub.model]) == "NUCLEOTIDE"){
    ret <- matrix(ret, nrow = 4, ncol = 4, byrow = TRUE)
  } else{
    ret <- matrix(ret, nrow = 2, ncol = 2, byrow = TRUE)
  }
  attr(ret, "class") <- "Pt"
  attr(ret, "code.type") <- as.character(.substitution$code.type[sub.model])
  attr(ret, "log") <- log
  ret
}

### Print Pt matrix.
print.Pt <- function(x, ...){
  Pt <- x
  if(attr(Pt, "code.type") == "NUCLEOTIDE"){
    colnames(Pt) <- .nucleotide$code[1:4]
    rownames(Pt) <- .nucleotide$code[1:4]
  } else{
    colnames(Pt) <- .snp$code[1:2]
    rownames(Pt) <- .snp$code[1:2]
  }
  if(attr(Pt, "log") == TRUE){
    cat("logPt\n")
  } else{
    cat("Pt\n")
  }
  attr(Pt, "class") <- NULL
  attr(Pt, "code.type") <- NULL
  attr(Pt, "log") <- NULL
  my.print(Pt)
}

