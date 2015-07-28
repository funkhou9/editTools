#' Processes VCF files in search for candidate RNA editing events
#' 
#' Must supply two files - one from RNA seq alignments that come from
#'  plus strand transcripts (VCF file 1). The other from minus strand transcripts.
#'  (VCF file 2).
#'  
#' Each vcf file requires two samples, a genomic DNA sample (from a WGS bam file)
#'  and a RNA sample coming from a tissue of interest (from a RNAseq bam file)
#'
#' @param file_plus input filename for VCF file 1.
#' @param file_minus input filename for VCF file 2.
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
                       file_minus,
                       qual = 10,
                       ex_indel = TRUE,
                       geno_dp = 10,
                       geno_hom = 95,
                       edit_dp = 5) {
  
  # Initialize lh references to chars
  p <- '+'
  m <- '-'
  
  # Process files by - 
  # 1. capture stdout,
  # 2. split by tab delimiter,
  # 3. store into matrix
  # 4. conver to df
  
  # Process plus file
  plus <- capture.output(edit_search(file_plus,
                                     p,
                                     qual,
                                     ex_indel,
                                     geno_dp,
                                     geno_hom,
                                     edit_dp)) %>%
            sapply(function(x) strsplit(x, split = '\t')) %>%
              do.call(rbind, .) %>%
                as.data.frame(stringsAsFactors = FALSE)
  
              
  
  # Process minus file
  minus <- capture.output(edit_search(file_minus,
                                      m,
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
  
  return (result)
}