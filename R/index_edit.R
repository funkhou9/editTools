# Indexes edits
#
# @param snp A single $sample entry
index_edit <- function(snp) {
  dna <- as.character(snp[1])
  rna <- as.character(snp[2])
  
  if(dna == "0/0") dna <- 1
  if(dna == "1/1") dna <- 2
  if(rna == "0/0") rna <- 1
  if(rna == "1/1") rna <- 2
  
  if(rna == "0/1") {
    if(dna == 1) rna <- 2
  } else if(dna == 2) rna <- 1
  
  return(data.frame(dna, rna))
}