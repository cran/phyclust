### This file contains functions for "ms".

ms <- function(nsam = NULL, nreps = 1, opts = NULL){
  temp.file <- tempfile("ms.")

  if(! is.null(opts) && ! is.null(nsam)){
    if(nsam >= 2){
      nsam <- as.character(nsam)
      nreps <- as.character(nreps)
      argv <- c("ms", nsam, nreps, unlist(strsplit(opts, " ")))

      .Call("R_ms_main", argv, temp.file, PACKAGE = "phyclust")
      ret <- scan(file = temp.file, what = "character", sep = "\n", quiet = TRUE)
      unlink(temp.file)

      class(ret) <- "ms"
      return(ret)
    }
  }

  argv <- c("ms", "-h")
  try(.Call("R_ms_main", argv, temp.file, PACKAGE = "phyclust"), silent = TRUE)
  unlink(temp.file)
  invisible()
} # End of ms().

print.ms <- function(x, ...){
  ms <- x
  cat(ms, sep = "\n")
} # End of print.ms().
