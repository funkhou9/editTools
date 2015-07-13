# Helpful accessor functions 
#
# @param vcf objects
# @return either SNPs or Samples component of vcf object
snps <- function(obj) obj$SNPs
samples <- function(obj, n)  obj$Samples[n]


#' Subsetting method for vcf class
#' 
#' @param obj a vcf obect
#' @param p numeric to filter SNPs
#' @param n numeric to filter Samples
#' @return a subsetted vcf object
#' @export
`[.vcf` <- function(obj, p, n) {
  #Method for `[` that will be used in other functions to filter and subset
  L <- list(SNPs = snps(obj)[p,], 
            Samples = lapply(samples(obj, n), function(x){x[p, ]}))
  class(L) <- "vcf"
  return(L)
}


# Merging tool to bring together genotypes (homo var, het, or homo ref) from
# different samples together.
#
# @param object a vcf object
# @param number of samples in vcf file
# @param col which columns would you like to merge from each sample? examples:1,2,3 or "GT", "PL"
# @return a merged sample object
SampleMerge <- function(obj, NS, col) {
  cols <- lapply(samples(obj, 1:NS), function(x){data.frame(x[,col])})
  call <- do.call("cbind", cols)
  return(call)
}


# Indexes genotypes
#
# @param x a single genotype - "0/0", "0/1" or "1/1"
indexGeno <- function(x) {
  if (x == "0/0") return(1)
  if (x == "1/1") return(2)
  if (x == "0/1") return(2) else return(NA)
}


# Indexes edits
#
# @param snp A single $sample entry
indexEdit <- function(snp) {
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


# Take in a vcf class and return genotype calls for each sample.
#
# creates an error when, after filters, there is nothing left! Rather than reporting
# an error, would like to report a message like "no edits found" 
#
# @param obj vcf object
# @param number of samples in vcf file (required)
# @return Mismatch results!
callGeno <- function(obj, NS) {
  
  # Merging samples to compare
  calls <- SampleMerge(obj, NS, 1)
  # Extract reference and alternative calls
  geno <- snps(obj)[, c("REF", "ALT")]
  
  # Index edits. Use to prepare a matrix for use in indexing
  i <- apply(calls, 1, indexEdit)
  idf <- data.frame(matrix(unlist(i), ncol = NS, byrow = TRUE))
  mat <- lapply(idf, function(x) cbind(seq_along(x), x))
  genocalls <- lapply(mat, function(x){geno[x]})
  genocalls <- data.frame(matrix(unname(unlist(genocalls)), ncol=NS))
  
  names(genocalls) <- names(samples(obj))
  genocalls <- cbind(snps(obj)[, c("CHROM", "POS")], genocalls, SampleMerge(obj, NS, 1:4))
  genocallcounts <- data.frame(table(genocalls[, names(samples(obj))]))
  freq <- ncol(genocallcounts)
  genocallcounts <- genocallcounts[genocallcounts[, freq] != 0, ]
  tot <- sum(genocallcounts[, freq])
  prob <- round(as.numeric(as.character(genocallcounts[, freq]))/tot, 3) 
  genocallcounts <- cbind(genocallcounts, prob)
  genocallcounts <- genocallcounts[with(genocallcounts, order(Freq, decreasing=T)), ]
  return(list(Calls = genocalls, Summary = genocallcounts))
}