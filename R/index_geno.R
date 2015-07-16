# Indexes genotypes
#
# @param x a single genotype - "0/0", "0/1" or "1/1"
# @return index {1, 2, NA} corresponding to genotype
index_geno <- function(x) {
#   if (x == "0/0") return(1)
#   if (x == "1/1") return(2)
#   if (x == "0/1") return(2) else return(NA)
  x <- switch(x,
              "0/0" = 1,
              "0/1" = 2,
              "1/1" = 2)

}