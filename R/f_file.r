### Functions to read and write files.

read.fasta <- function(filename, byrow = TRUE, code.type = .code.type[1]){
  if(code.type[1] == "NUCLEOTIDE"){
    ret <- read.fasta.nucleotide(filename, byrow = byrow)
  } else if(code.type[1] == "SNP"){
    stop("SNP is not implemented yet.")
  } else{
    stop("code.type is not found.")
  }

  ret
} # End of read.fasta().

read.phylip <- function(filename, byrow = TRUE, code.type = .code.type[1]){
  if(code.type[1] == "NUCLEOTIDE"){
    ret <- read.phylip.nucleotide(filename, byrow = byrow)
  } else if(code.type[1] == "SNP"){
    ret <- read.phylip.snp(filename, byrow = byrow)
  } else{
    stop("code.type is not found.")
  }

  ret
} # End of read.phylip().

write.fasta <- function(seqdata, filename, classid = NULL,
    seqname = NULL, width.line = 60, lower.case = FALSE,
    code.type = .code.type[1]){
  if(code.type[1] == "NUCLEOTIDE"){
    ret <- write.fasta.nucleotide(seqdata, filename, classid = classid,
             seqname = seqname, width.line = width.line,
             lower.case = lower.case)
  } else if(code.type[1] == "SNP"){
    stop("SNP is not implemented yet.")
  } else{
    stop("code.type is not found.")
  }
} # End of write.fasta().

write.phylip <- function(seqdata, filename, classid = NULL,
    seqname = NULL, width.seqname = 10, width.line = 60, lower.case = FALSE,
    code.type = .code.type[1]){
  if(code.type[1] == "NUCLEOTIDE"){
    ret <- write.phylip.nucleotide(seqdata, filename, classid = classid,
             seqname = seqname, width.seqname = width.seqname,
             width.line = width.line, lower.case = lower.case)
  } else if(code.type[1] == "SNP"){
    ret <- write.phylip.snp(seqdata, filename, classid = classid,
             seqname = seqname, width.seqname = width.seqname,
             width.line = width.line)
  } else{
    stop("code.type is not found.")
  }
} # End of write.phylip().

