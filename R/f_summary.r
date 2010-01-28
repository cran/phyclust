### This file contains functions for print and summary.

### Print phyclust
print.phyclust <- function(x, digits = max(4, getOption("digits") - 3), ...){
  pcobj <- x

  options(digits = digits)
  init <- NULL

  if(!is.null(pcobj$conv)){
    my.cat("Phyclust Results:\n",
           "code type: ", pcobj$code.type,
           ", em method: ", pcobj$em.method,
           ", boundary method: ", pcobj$boundary.metho, ".\n",
           "init procedure: ", pcobj$init.procedure,
           ", method: ", pcobj$init.method, ".\n",
           "model substitution: ", pcobj$substitution.model,
           ", distance: ", pcobj$edist.model, ".\n",
           "iter: ", pcobj$conv$iter,
           " ", pcobj$conv$inner.iter,
           " ", pcobj$conv$cm.iter,
           ", convergence: ", pcobj$conv$flag,
           ", check.param: ", pcobj$conv$check.param, ".\n",
           "eps: ", pcobj$conv$eps,
           ", error: ", pcobj$conv$error, ".\n")
  }
      
  my.cat("N.X.org: ", pcobj$N.X.org,
         ", N.X.unique: ", pcobj$N.X.unique,
         ", L: ", pcobj$L,
         ", K: ", pcobj$K,
         ", p: ", pcobj$p,
         ", N.seg.site: ", pcobj$N.seg.site, ".\n",
         "logL: ", pcobj$logL)
  if(!is.null(pcobj$bic)) my.cat(", bic: ", pcobj$bic)
  if(!is.null(pcobj$aic)) my.cat(", aic: ", pcobj$aic)
  if(!is.null(pcobj$icl)) my.cat(", icl: ", pcobj$icl)
  my.cat("\n")
  my.cat("identifier: ", pcobj$QA$identifier, "\n")
  cat("  Eta:", pcobj$Eta, "\n")
  if(!is.null(pcobj$Q$pi)){
    my.cat("  pi:\n")
    my.print(pcobj$Q$pi)
  }
  if(!is.null(pcobj$Q$kappa)) cat("  kappa:", pcobj$Q$kappa, "\n")
  cat("  Tt:", pcobj$Q$Tt, "\n")
  if(!is.null(pcobj$n.class)) cat("  n.class:", pcobj$n.class, "\n")
} # End of print.phyclust().

summary.phyclust <- function(object, ...){
  pcobj <- object

  print.phyclust(pcobj)
  cat("Mu:\n")
  for(i in 1:pcobj$K){
    cat("   ", as.character(.nucleotide$code[pcobj$Mu[i,] + 1]), "\n")
  }
  cat("Class id:", pcobj$class.id, "\n")
} # End of summary.phyclust().

