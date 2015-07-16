# Take in a vcf class and return genotype calls for each sample.
#
# creates an error when, after filters, there is nothing left! Rather than reporting
# an error, would like to report a message like "no edits found" 
#
# @param obj vcf object
# @param number of samples in vcf file (required)
# @return Mismatch results!
#' @export
call_geno <- function(obj, NS) {
  
  # Merging sample calls to compare
  calls <- sample_merge(obj, NS, 1)
  
  # Extract reference and alternative calls
  geno <- snps(obj)[, c("REF", "ALT")]

  # Index edits. Use to prepare a matrix for use in indexing
  i <- apply(calls, 1, index_edit)
  
  idf <- data.frame(matrix(unlist(i), ncol = NS, byrow = TRUE))
  
  
  mat <- lapply(idf, function(x) cbind(seq_along(x), x))
  genocalls <- lapply(mat, function(x){geno[x]})
  genocalls <- data.frame(matrix(unname(unlist(genocalls)), ncol=NS))
  
  
  
  names(genocalls) <- names(samples(obj))
  genocalls <- cbind(snps(obj)[, c("CHROM", "POS")], genocalls, sample_merge(obj, NS, 1:4))
  genocallcounts <- data.frame(table(genocalls[, names(samples(obj))]))
  freq <- ncol(genocallcounts)
  genocallcounts <- genocallcounts[genocallcounts[, freq] != 0, ]
  tot <- sum(genocallcounts[, freq])
  prob <- round(as.numeric(as.character(genocallcounts[, freq]))/tot, 3) 
  genocallcounts <- cbind(genocallcounts, prob)
  genocallcounts <- genocallcounts[with(genocallcounts, order(Freq, decreasing=T)), ]
  return(list(Calls = genocalls, Summary = genocallcounts))
}