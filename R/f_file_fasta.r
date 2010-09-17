### Read FASTA file.
### Return a list contains
###   nseq: Number of sequence.
###   seqlen: Sequence length.
###   seqname: Sequence's name.
###   org: Original sequence in nid. array[nsite, nseq]

read.fasta.nucleotide <- function(filename, byrow = TRUE){
  tmp <- readLines(filename)
  org <- NULL
  tmp.org <- NULL
  seqname <- NULL

  i <- 1
  nseq <- 0
  repeat{
    if(regexpr("^>", tmp[i]) == 1){
      if(!is.null(tmp.org)){
        org <- cbind(org, unlist(strsplit(tmp.org, "")))
        tmp.org <- NULL
      }

      if(regexpr("\\|", tmp[i]) == 1){
        tmp.seqname <- unlist(strsplit(tmp[i], "\\|"))[1]
      } else{
        tmp.seqname <- unlist(strsplit(tmp[i], " "))[1]
      }

      seqname <- c(seqname, gsub(">(.*)", "\\1", tmp.seqname))
      nseq <- nseq + 1
    } else if(tmp[i] %in% c("", " ")){
    } else{
      tmp.org <- paste(tmp.org, tmp[i], sep = "")
    }

    i <- i + 1
    if(i > length(tmp)){
      org <- cbind(org, unlist(strsplit(tmp.org, "")))
      break
    }
  }
  
  fasta <- list(code.type = "NUCLEOTIDE",
              nseq = nseq, seqlen = nrow(org), seqname = seqname,
              org.code = org, org = NULL, byrow = byrow)
  fasta$org <- matrix(code2nid(org), nrow = fasta$seqlen, ncol = fasta$nseq)

  if(byrow){
    fasta$org.code <- t(fasta$org.code)
    fasta$org <- t(fasta$org)
  }

  class(fasta) <- "seq.data"
  fasta
} # End of read.fasta.nucleotide().


write.fasta.nucleotide <- function(seqdata, filename, classid = NULL,
    seqname = NULL, width.line = 60, lower.case = FALSE){
  seqdata <- apply(seqdata, 2, nid2code, lower.case = lower.case)
  
  n.seq <- nrow(seqdata)
  tl.seq <- ncol(seqdata)

  if(is.null(seqname)){
    seqname <- as.character(1:n.seq)
  }

  if(!is.null(classid)){
    seqname <- cbind(seqname, rep("-", n.seq) , classid)
    seqname <- apply(seqname, 1, paste, collapse = "")
  }

  tl.show <- ceiling(tl.seq / width.line)
  show.range <- list()
  for(i in 1:tl.show){
    show.range[[i]] <- (1:width.line) + (i - 1) * width.line
    if(show.range[[i]][width.line] > tl.seq){
      show.range[[i]] <- show.range[[i]][1:which(show.range[[i]] == tl.seq)]
    }
  }

  ret <- NULL
  for(i in 1:n.seq){
    ret <- c(ret, paste(">", seqname[i], sep = ""))
    for(j in 1:tl.show){
      ret <- c(ret, paste(seqdata[i, show.range[[j]]], collapse = ""))
    }
  }

  write.table(ret, file = filename, quote = FALSE, sep = "", row.names = FALSE,
              col.names = FALSE)
} # End of write.fasta.nucleotide().

