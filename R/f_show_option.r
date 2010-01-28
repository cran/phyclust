### This file contains function to show all possible options in phyclust.

.show.option <- function(){
  my.cat("Boundary method: ", paste(.boundary.method, collapse = ", "), "\n")
  my.cat("Code type: ", paste(.code.type, collapse = ", "), "\n")
  my.cat("Edist model: ", paste(.edist.model, collapse = ", "), "\n")
  my.cat("EM method: ", paste(.em.method, collapse = ", "), "\n")
  my.cat("Identifier: ", paste(.identifier, collapse = ", "), "\n")
  my.cat("Init method: ", paste(.init.method, collapse = ", "), "\n")
  my.cat("Init procedure: ", paste(.init.procedure, collapse = ", "), "\n")
  my.cat("Standard code: \n")
  my.print(as.matrix(.nucleotide))
  my.print(as.matrix(.snp))
  my.cat("Substitution model: \n")
  my.print(as.matrix(.substitution))
  invisible()
} # End of .show.option()
