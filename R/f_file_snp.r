### Read SNP file. (PHYLIP format)
### Return a list contains
###   nseq: Number of sequence.
###   seqlen: Sequence length.
###   seqname: Sequence's name.
###   org: Original sequence in nid. array[nseq, nsite]
read.phylip.snp <- function(filename, byrow = TRUE){
  snp <- list(code.type = "SNP", info = NULL, nseq = NULL,
              seqlen = NULL, seqname = NULL, org.code = NULL, org = NULL,
              byrow = byrow)

### Read header.
  snp$info <- readLines(filename, n = 1)
  tmp <- unstrsplit(snp$info, " ")
  tmp <- tmp[tmp != ""]
  snp$nseq <- as.numeric(tmp[1])
  snp$seqlen <- as.numeric(tmp[2])

### Read data and transfer to nid.
  op <- options("stringsAsFactors")
  options(stringsAsFactors = FALSE)
  tmp <- read.table(filename, sep = "", quote = "", skip = 1, fill = TRUE)
  options(op)

### Split the data by reading blocks and rejoin them by sequences.
  org <- split(tmp, gl(nrow(tmp) / snp$nseq, snp$nseq))
  org <- do.call("cbind", org)
  org <- as.matrix(org, nrow = snp$nseq)
  snp$seqname <- org[, 1]
  org <- apply(as.matrix(org[, -1], nrow = nrow(org)), 1,
               function(y){ unstrsplit(as.character(y), "") })
  snp$org.code <- matrix(org, nrow = snp$seqlen, ncol = snp$nseq)
  snp$org <- matrix(snp2sid(org), nrow = snp$seqlen, ncol = snp$nseq)

  if(byrow){
    snp$org.code <- t(snp$org.code)
    snp$org <- t(snp$org)
  }

  class(snp) <- "seq.data"
  snp
} # End of read.phylip.snp().


write.phylip.snp <- function(seqdata, filename, classid = NULL,
    seqname = NULL, width.seqname = 10, width.line = 60){
  seqdata <- apply(seqdata, 2, sid2snp)

  n.seq <- nrow(seqdata)
  tl.seq <- ncol(seqdata)

  if(is.null(seqname)){
    seqname <- as.character(1:n.seq)
  }

  if(!is.null(classid)){
    seqname <- cbind(seqname, rep("-", n.seq) , classid)
    seqname <- apply(seqname, 1, paste, collapse = "")
  }

  tmp.seqname <- strsplit(seqname, "")
  for(i in 1:n.seq){
    tl <- length(tmp.seqname[[i]])
    if(tl > width.seqname){
      tmp.seqname[[i]] <- paste(tmp.seqname[[i]][1:width.seqname],
                                collapse = "")
    } else{
      tmp.seqname[[i]] <- paste(c(tmp.seqname[[i]],
                                  rep(" ", width.seqname - tl)), collapse = "")
    }
  }
  tmp.seqname <- do.call("c", tmp.seqname)

  head <- paste(n.seq, tl.seq, collapse = " ")
  write(head, file = filename)

  tl.show <- ceiling(tl.seq / width.line)
  for(i in 1:tl.show){
    show.range <- (1:width.line) + (i - 1) * width.line
    if(show.range[width.line] > tl.seq){
      show.range <- show.range[1:which(show.range == tl.seq)]
    }

    tmp.seqdata <- apply(seqdata[, show.range], 1, paste, collapse = "")

    if(i == 1){
      tmp.seqdata <- cbind(tmp.seqname, tmp.seqdata)
      tmp.seqdata <- apply(tmp.seqdata, 1, paste, collapse = "")
    } 

    write.table(tmp.seqdata, file = filename, append = TRUE, col.names = FALSE,
                row.names = FALSE, sep = "", quote = FALSE)
    write("", file = filename, append = TRUE)
  }
} # End of write.phylip.snp().

