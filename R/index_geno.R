# @param x a vector of genotypes - assumes x {"0/0", "0/1", "1/1"}
# @return a vector of genotype indexes {1, 2}
index_geno <- function(x) {
  
  # Homozygous reference (0/0) get 1 coding
  #   Others  (0/1, 1/1) get 2 coding
  idx_geno <- c(1, 2)[(x == "1/1") + 1]
  
  return (idx_geno)
}