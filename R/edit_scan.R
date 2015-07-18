#' Performs a scan of all putative RNA editing sites from vcf files
#' 
#' Requires two VCF files as input. One of the files requires RNA sequencing data
#'  from PLUS strand alignments, the other from MINUS strand alignments.
#'  Each file needs to have two samples. One sample will always be from genomic DNA. 
#'  The other will be from a RNA sequencing data from a tissue of interest.
#'  
#' @param plus_file
#' @param minus_file
#' @param qual An integer specifiying the minimum variant QUAL
#' @param ex.indel logical indicating whether to exclude indels from the scan
#' @param geno.dp integer specifying the minimum genotype depth
#' @param geno.hom integer ranging from 0 to 1 specifiying the proportion of homozygosity
#'  the genotype must exhibit
#' @param edit.dp integer specifying the minimum depth required for evidence of
#'  RNA editing
#' @param summary logical indicating whether a summary of found edits is returned
#'  or a subsetted vcf object
#' @return an edit_summary object
#' @export
edit_scan <- function(plus_file,
                      minus_file,
                      names = c("DNA", "RNA"),
                      qual = 10,
                      ex.indel = TRUE,
                      geno.dp = 10,
                      geno.hom = 95,
                      edit.dp = 5,
                      summary = TRUE) {
  plus_strand <- 
    read_vcf(plus_file, names) %>%
      edit_summary(qual_ = qual,
                   ex.indel_ = ex.indel,
                   geno.dp_ = geno.dp,
                   geno.hom_ = geno.hom,
                   edit.dp_ = edit.dp,
                   summary_ = summary)
  
  minus_strand <-
    read_vcf(minus_file, names) %>%
    edit_summary(qual_ = qual,
                 ex.indel_ = ex.indel,
                 geno.dp_ = geno.dp,
                 geno.hom_ = geno.hom,
                 edit.dp_ = edit.dp,
                 summary_ = summary)
  
  merged_results <- merge_strands(plus_strand, minus_strand)
} 