### Object Name (Type): Definaiton (Examples)
###   codeseq (string): A nucleotide sequence. (ACGT...)
###   nidseq (vector[L]): A nucleotide id sequence. (1234....)
###                       The length should be L = 2 + 3 * n + 3 for n >= 1.


### Transfer an nucleotide sequence to a code id sequence.
### Eg. ACGT-... => 01234...
code2nid <- function(codeseq){
  nidseq <- codeseq
  for(i in 1:nrow(.nucleotide)){
    nidseq[codeseq == as.character(.nucleotide$code[i])] <- .nucleotide$nid[i]
    nidseq[codeseq == as.character(.nucleotide$code.l[i])] <- .nucleotide$nid[i]
  }
  nidseq[! codeseq %in%
         c(as.character(.nucleotide$code), as.character(.nucleotide$code.l))] <- 4
  nidseq <- as.numeric(nidseq)

  if(is.matrix(codeseq)){
    dim(nidseq) <- dim(codeseq)
  }
  nidseq
} # End of code2nid().

### Transfer a nucleotide id sequence to an nucleotide sequence.
### Eg. 01234... => ACGT-...
nid2code <- function(nidseq, lower.case = TRUE){
  codeseq <- nidseq
  for(i in 1:nrow(.nucleotide)){
    if(lower.case){
      codeseq[nidseq == .nucleotide$nid[i]] <- as.character(.nucleotide$code.l[i])
    } else{
      codeseq[nidseq == .nucleotide$nid[i]] <- as.character(.nucleotide$code[i])
    }
  }
  codeseq[!nidseq %in% .nucleotide$nid] <- "-"

  if(is.matrix(nidseq)){
    dim(codeseq) <- dim(nidseq)
  }
  codeseq
} # End of nid2code().

### Transfer an SNP sequence to a SNP id sequence.
### Eg. 1212-... => 01012...
snp2sid <- function(snpseq){
  sidseq <- snpseq
  for(i in 1:nrow(.snp)){
    sidseq[snpseq == as.character(.snp$code[i])] <- .snp$sid[i]
  }
  sidseq[! snpseq %in% as.character(.snp$code)] <- 2
  sidseq <- as.numeric(sidseq)

  if(is.matrix(snpseq)){
    dim(sidseq) <- dim(snpseq)
  }
  sidseq
} # End of snp2sid().

### Transfer a SNP id sequence to an SNP sequence.
### Eg. 01012... => 1212-..
sid2snp <- function(sidseq){
  snpseq <- sidseq
  for(i in 1:nrow(.snp)){
    snpseq[sidseq == .snp$sid[i]] <- as.character(.snp$code[i])
  }
  snpseq[!sidseq %in% .snp$sid] <- "-"

  if(is.matrix(sidseq)){
    dim(snpseq) <- dim(sidseq)
  }
  snpseq
} # End of sid2snp().


### Transfer a nucleotide seqeunce to a SNP sequence.
### Eg. AGCT-... => 1122-...
code2snp <- function(codeseq){
  snpseq <- codeseq
  snpseq[codeseq %in% c("A", "a", "G", "g")] <- "1"
  snpseq[codeseq %in% c("C", "c", "T", "t")] <- "2"
  snpseq[!(snpseq %in% c("1", "2"))] <- "-"

  if(is.matrix(codeseq)){
    dim(snpseq) <- dim(codeseq)
  }
  snpseq
} # End of code2snp().

### Transfer a SNP seqeunce to a nucleotide sequence.
### Eg. 1122-... => AACC-...
snp2code <- function(snpseq, half = TRUE){
  codeseq <- snpseq

  if(half){
    id <- snpseq == "1"
    tl <- sum(id)
    codeseq[id] <- rep(c("A", "G"), ceiling(tl / 2))[1:tl]
    id <- snpseq == "2"
    tl <- sum(id)
    codeseq[id] <- rep(c("C", "T"), ceiling(tl / 2))[1:tl]
  } else{
    codeseq[snpseq == "1"] <- "A"
    codeseq[snpseq == "2"] <- "C"
  }

  codeseq[!(snpseq %in% c("1", "2"))] <- "-"

  if(is.matrix(snpseq)){
    dim(codeseq) <- dim(snpseq)
  }
  codeseq
} # End of code2snp().


### Transfer nid and sid.
nid2sid <- function(nidseq){
  snp2sid(code2snp(nid2code(nidseq)))
} # End of nid2sid()

sid2nid <- function(sidseq, half = TRUE){
  code2nid(snp2code(sid2snp(sidseq), half = half))
} # End of nid2sid()
