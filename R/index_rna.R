# @param rna a vector of RNA calls {"0/0", "0/1", "1/1"}
# @param index_geno a vector of genotype indexes produced by index_geno
# @return a vector of RNA indexes {1, 2}
index_rna <- function(rna, index_geno) {
  
  # If RNA is heterozygous - need to look at DNA
  if (rna == "0/1") {
    idx_rna <- c(1, 2)[(index_geno == 1) + 1]
  } else {
    idx_rna <- c(1, 2)[(rna == "1/1") + 1]
  }
}




