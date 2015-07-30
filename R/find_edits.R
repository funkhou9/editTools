#' Processes VCF files in search for candidate RNA editing events
#' 
#' Must supply two files - one from RNA seq alignments that come from
#'  plus strand transcripts (VCF file 1). The other from minus strand transcripts.
#'  (VCF file 2).
#'  
#' Each vcf file requires a genomic DNA sample (the first sample listed in the vcf file),
#'  along with any number of RNA samples from various tissues.
#'
#' @param file_plus input filename for VCF file 1.
#' @param file_minus input filename for VCF file 2. If only one VCF file provided, it is assumed
#'  that the RNA variants in that file all come from plus strand transcripts.
#' @param names A character vector 
#' @param qual An integer specifiying the minimum variant QUAL
#' @param ex_indel logical indicating whether to exclude indels from the scan
#' @param geno_dp integer specifying the minimum genotype depth
#' @param geno_hom integer ranging from 0 to 1 specifiying the proportion of homozygosity
#'  the genotype must exhibit
#' @param edit_dp integer specifying the minimum depth required for evidence of
#'  RNA editing
#' @param summary logical indicating whether a summary of found edits is returned
#'  or a subsetted vcf object
#' @return an edit_summary object
#' @import magrittr
#' @export
find_edits <- function(file_plus,
                       file_minus = NULL,
                       names = character(),
                       qual = 10,
                       ex_indel = TRUE,
                       geno_dp = 10,
                       geno_hom = 95,
                       edit_dp = 5) {
  
  # Initialize lh references to chars
  p <- "+"
  m <- "-"
  
  # Requre tissue samples to be labeled by user
  if (length(names) == 0) {
    stop ("Please provide names argument")
  }
  
  # Process files by - 
  # 1. capture stdout,
  # 2. split by tab delimiter,
  # 3. store into matrix
  # 4. conver to df
  
  # Process plus file
  plus <- capture.output(edit_search(file_plus,
                                     p,
                                     names,
                                     qual,
                                     ex_indel,
                                     geno_dp,
                                     geno_hom,
                                     edit_dp)) %>%
            sapply(function(x) strsplit(x, split = '\t')) %>%
              do.call(rbind, .) %>%
                as.data.frame(stringsAsFactors = FALSE)
  
              
  if (!is.null(file_minus)) {
    # Process minus file
    minus <- capture.output(edit_search(file_minus,
                                        m,
                                        names,
                                        qual,
                                        ex_indel,
                                        geno_dp,
                                        geno_hom,
                                        edit_dp)) %>%
      sapply(function(x) strsplit(x, split = '\t')) %>%
      do.call(rbind, .) %>%
      as.data.frame(stringsAsFactors = FALSE)
    
    # Some df formatting -
    # 1. Combine plus strand and minus strand results
    # 2. Format sample_names arg into a usable header, apply header
    # 3. Enusre POS is numeric, then reorder by CHROM, then POS
    # 4. Remove rownames
    result <- rbind(plus, minus)
  } else
    result <- plus

  colnames(result) <- c("CHR",
                        "POS",
                        "Strand",
                        "DNA",
                        "RNA",
                        "DNA_DP",
                        "DNA_DV",
                        "RNA_DP",
                        "RNA_DV",
                        "Tissue")
  
  result[, "POS"] <- as.numeric(result[, "POS"])
  result <- result[order(result[, "CHR"], result[,"POS"]), ]

  rownames(result) <- NULL
  
  # Produce mismatch counts for each tissue sample
  mismatches <- lapply(names,
                       count_mismatch,
                       result)
  names(mismatches) <- names
  
  # Append mismatch counts to existing results and declare class
  result <- list("Edits" = result,
                 "Mismatches" = mismatches)
  
  class(result) <- "edit_table"
  return (result)
}