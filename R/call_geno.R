# Take in a vcf class and return genotype calls for each sample.
#
# creates an error when, after filters, there is nothing left! Rather than reporting
# an error, would like to report a message like "no edits found" 
#
# @param obj vcf object
# @param number of samples in vcf file (required)
# @return Mismatch results!
# @export
call_geno <- function(obj, NS) {
  
  # Merging sample calls to compare
  calls <- sample_merge(obj, NS, 1)
  
  # Extract reference and alternative calls
  geno <- snps(obj)[, c("REF", "ALT")]

  # Index DNA genotypes and prepare for use in matrix indexing
  idx_dna <- index_geno(as.character(calls[, 1]))
  
  # For each locus at a time, evaluate rna index according to
  #   dna index
  # MIGHT BE ABLE TO CHANGE THIS STEP TO ACCOMODATE MORE SAMPLES
  idx_rna <- mapply(index_rna,
                    as.character(calls[, 2]),
                    idx_dna)
  
  # Prepare each index to be used as an index
  idx_dna <- cbind(seq_along(idx_dna), idx_dna)
  idx_rna <- cbind(seq_along(idx_rna), idx_rna)
  
  # Obtain DNA and RNA calls using indexes and merge
  dna_calls <- geno[idx_dna]
  rna_calls <- geno[idx_rna]
  
  dna_rna_calls <- cbind(dna_calls, rna_calls)
  
  # Give results names of samples - i.e. "DNA", "RNA"
  colnames(dna_rna_calls) <- names(samples(obj))

  # Reformat to include additional information from "snp" and
  #   "sample" fields
  dna_rna_calls <- cbind(snps(obj)[, c("CHROM", "POS")],
                         dna_rna_calls,
                         sample_merge(obj, NS, 1:4))
  
  # Obtain counts of each DNA/RNA mismatch
  edit_counts <-
    table(dna_rna_calls[, names(samples(obj))]) %>%
      data.frame()
  
  # Remove "Zero" mismatch types
  edit_counts <- edit_counts[edit_counts[, "Freq"] != 0, ]
  
  # Get total mismatches
  mismatch_tot <- sum(edit_counts[, "Freq"])
  
  # Proportions of total mismatches
  mismatch_prob <- 
    (edit_counts[, "Freq"] / mismatch_tot) %>%
      round(3)
  edit_counts <- cbind(edit_counts, "Prob" = mismatch_prob)
  
  # Reorder so that most common mismatch is on top
  edit_counts <- edit_counts[with(edit_counts, order(Freq, decreasing=T)), ]
  
  result <- list("Calls" = dna_rna_calls,
                 "Summary" = edit_counts)
  
  class(result) <- "edit_summary"
  return(result)
}