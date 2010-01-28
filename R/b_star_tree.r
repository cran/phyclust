### This file contains functions for PARAMETERIC BOOTSTRAP that generate
### star trees based on an pcobj fitted by phyclust(), and generate sequences
### according to the star trees.

bootstrap.star.tree <- function(n.tip, height = NULL){
  if(n.tip > 1){
    ms.tree <- ms(n.tip, opts = "-T")
    rooted.tree <- read.tree(text = ms.tree[3])
    star.tree <- as.star.tree(rooted.tree)
  } else if(n.tip == 1){
    star.tree <- list(edge = matrix(c(3, 3, 1, 2), nrow = 2),
                      Nnode = 1,
                      tip.label = c("1", "2"),
                      edge.length = c(1, 0))
    class(star.tree) <- "phylo"
  } else{
    stop("n.tip > 0")
  }

  if(! is.null(height)){
    star.tree$edge.length[star.tree$edge[, 2] <= n.tip] <- height
  }
  star.tree$n.tip <- n.tip

  star.tree
} # End of bootstrap.star.tree().

bootstrap.star.trees <- function(pcobj, min.n.class = 1){
  if(min.n.class * pcobj$K >= pcobj$N.X.org){
    stop("min.n.class * pcobj$k >= pcobj$N.X.org")
  }

  repeat{
    n.class <- rmultinom(1, pcobj$N.X.org, pcobj$Eta)
    if(all(n.class >= min.n.class)){
      n.class <- as.vector(n.class)
      break
    }
  }

  Tt <- pcobj$QA$Tt
  if(length(Tt) == 1){		# EE or VE
    Tt <- rep(Tt, pcobj$K)
  }

  star.trees <- NULL
  for(k in 1:pcobj$K){
    star.trees[[k]] <- bootstrap.star.tree(n.class[k], height = Tt[k])
  }

  star.trees
} # End of bootstrap.star.trees().

bootstrap.merge.seq <- function(seqs, code.type = .code.type[1]){
  K <- length(seqs)
  da <- read.seqgen(seqs[[1]], code.type = code.type[1])
  if(K > 1){
    for(k in 2:K){
      tmp.da <- read.seqgen(seqs[[k]], code.type = code.type[1])
      da$nseq <- da$nseq + tmp.da$nseq
      da$seqname <- c(da$seqname, tmp.da$seqname)
      da$org.code <- rbind(da$org.code, tmp.da$org.code)
      da$org <- rbind(da$org, tmp.da$org)
    }
    da$info <- paste(" ", da$nseq, " ", da$seqlen, sep = "")
  }
  da
} # End of bootstrap.merge.seq().

bootstrap.star.trees.seq <- function(pcobj, min.n.class = 1){
  star.trees <- bootstrap.star.trees(pcobj, min.n.class = min.n.class)
  seq.boot <- bootstrap.seq(pcobj, star.trees)
  list(trees = star.trees, seq = seq.boot)
} # End of bootstrap.star.trees.seq().

