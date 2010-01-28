### This file contains functions for PARAMETERIC BOOTSTRAP that generate
### star trees based on an pcobj fitted by phyclust(), and generate sequences
### according to the star trees.

bootstrap.seq <- function(pcobj, star.trees){
  if(pcobj$code.type == "NUCLEOTIDE"){
    seq.boot <- bootstrap.seq.nucleotide(pcobj, star.trees)
  } else if(pcobj$code.type == "SNP"){
    seq.boot <- bootstrap.seq.snp(pcobj, star.trees)
  } else{
    stop("code.type is not found.")
  }
  seq.boot
} # End of bootstrap.seq()

bootstrap.seq.nucleotide <- function(pcobj, star.trees){
  if(pcobj$code.type != "NUCLEOTIDE"){
    stop("Only for NUCLEOTIDE data.")
  }

  if(is.null(pcobj$QA$kappa)){
    kappa <- 1
  } else{
    kappa <- pcobj$QA$kappa
  }
  if(length(kappa) == 1){		# EE or EV
    kappa <- rep(kappa, pcobj$K)
  }

  if(is.null(pcobj$QA$pi)){
    pi <- rep(0.25, 4)
  } else{
    pi <- pcobj$QA$pi
  }
  if(length(pi) == 4){			# EE or EV
    pi <- matrix(rep(pi, pcobj$K), nrow = pcobj$K, byrow = TRUE)
  }

  seq.nucleo <- NULL
  for(k in 1:pcobj$K){
    seq.nucleo[[k]] <- gen.seq.HKY(star.trees[[k]], pi[k,], kappa[k],
                                   pcobj$L, anc.seq = pcobj$Mu[k,])
    if(star.trees[[k]]$n.tip == 1){
      seq.nucleo[[k]] <- seq.nucleo[[k]][1:2]
      tmp <- unstrsplit(seq.nucleo[[k]][1], " ")
      tmp <- tmp[tmp != ""]
      seq.nucleo[[k]][1] <- paste(" ", as.numeric(tmp[1]) - 1,
                                  " ", tmp[-1], sep = "")
      class(seq.nucleo[[k]]) <- "seqgen"
    }
  }

  seq.nucleo
} # End of bootstrap.seq.nucleotide().

bootstrap.seq.snp <- function(pcobj, star.trees){
  if(pcobj$code.type != "SNP"){
    stop("Only for SNP data.")
  }

  if(is.null(pcobj$QA$pi)){
    pi <- rep(0.5, 2)
  } else{
    pi <- pcobj$QA$pi
  }
  if(length(pi) == 2){			# EE or EV
    pi <- matrix(rep(pi, pcobj$K), nrow = pcobj$K, byrow = TRUE)
  }

  seq.snp <- NULL
  for(k in 1:pcobj$K){
    seq.snp[[k]] <- gen.seq.SNP(star.trees[[k]], pi[k,],
                                pcobj$L, anc.seq = pcobj$Mu[k,])
    if(star.trees[[k]]$n.tip == 1){
      seq.snp[[k]] <- seq.snp[[k]][1:2]
      tmp <- unstrsplit(seq.snp[[k]][1], " ")
      tmp <- tmp[tmp != ""]
      seq.snp[[k]][1] <- paste(" ", as.numeric(tmp[1]) - 1,
                               " ", tmp[-1], sep = "")
      class(seq.snp[[k]]) <- "seqgen"
    }
  }

  seq.snp
} # End of bootstrap.seq.snp().

