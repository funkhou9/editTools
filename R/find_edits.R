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
#' @param names A character vector specifying the names of RNA samples in the order they appear in the VCF file.
#' @param qual An integer specifiying the minimum variant QUAL
#' @param ex_indel logical indicating whether to exclude indels from the scan
#' @param geno_dp integer specifying the minimum genotype depth
#' @param geno_hom integer ranging from 0 to 1 specifiying the proportion of homozygosity
#'  the genotype must exhibit
#' @param edit_dp integer specifying the minimum depth required for evidence of
#'  RNA editing
#' @param strand_bias integer specifying maximum sample Phred-scaled strand bias for an RNA sample
#'  to be considered.
#' @return an edit_summary object
#' @import magrittr
#' @export
find_edits <- function(file_plus,
                       file_minus = NULL,
                       names = character(),
                       ex_indel = TRUE,
                       geno_dp = 10,
                       geno_hom = 95,
                       edit_dp = 5,
                       strand_bias = 20) {
  
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
  # 4. convert to df
  
  # Process plus file
  plus <- capture.output(edit_search(file_plus,
                                     p,
                                     names,
                                     ex_indel,
                                     geno_dp,
                                     geno_hom,
                                     edit_dp,
                                     strand_bias)) %>%
            sapply(function(x) strsplit(x, split = '\t')) %>%
              do.call(rbind, .) %>%
                as.data.frame(stringsAsFactors = FALSE)
  
              
  if (!is.null(file_minus)) {
    # Process minus file
    minus <- capture.output(edit_search(file_minus,
                                        m,
                                        names,
                                        ex_indel,
                                        geno_dp,
                                        geno_hom,
                                        edit_dp,
                                        strand_bias)) %>%
              sapply(function(x) strsplit(x, split = '\t')) %>%
                do.call(rbind, .) %>%
                  as.data.frame(stringsAsFactors = FALSE)
    
    # Some df formatting -
    # 1. Combine plus strand and minus strand results
    # 2. Ensure all numeric columns are numeric (Pos, DNA_depth:Ave_MQ), then reorder by CHROM, then POS
    # 3. Remove rownames
    result <- rbind(plus, minus)
  } else
    result <- plus

  # Add an "ID" column--doesn't do much. Just provides an identifier for a particular mismatch
  # found within a particular tissue. 
  result <- cbind(seq_along(result[, 1]), result)
  
  colnames(result) <- c("ID",
                        "Chr",
                        "Pos",
                        "Strand",
                        "Mismatch",
                        "DNA_depth",
                        "DNA_variant_depth",
                        "RNA_depth",
                        "RNA_mismatch_depth",
                        "RNA_edit_frac",
                        "Phred_strand_bias",
                        "Ave_MQ",
                        "Tissue")
  
  result[, "Pos"] <- as.numeric(result[, "Pos"])
  result[, "DNA_depth"] <- as.numeric(result[, "DNA_depth"])
  result[, "DNA_variant_depth"] <- as.numeric(result[, "DNA_variant_depth"])
  result[, "RNA_depth"] <- as.numeric(result[, "RNA_depth"])
  result[, "RNA_mismatch_depth"] <- as.numeric(result[, "RNA_mismatch_depth"])
  result[, "RNA_edit_frac"] <- as.numeric(result[, "RNA_edit_frac"])
  result[, "Phred_strand_bias"] <- as.numeric(result[, "Phred_strand_bias"])
  result[, "Ave_MQ"] <- as.numeric(result[, "Ave_MQ"])
  result <- result[order(result[, "Chr"], result[,"Pos"]), ]
  rownames(result) <- NULL
  
  # Convert columns intended to be numeric--
  
  result <- list("AllSites" = result)

  # Use count_mismatch() to get counts of each mismatch for each tissue
  mismatch_cnts <- count_mismatch(result$AllSites)
  
  # Append mismatch counts to existing results and declare class
  result <- append(result,
                   list("Tissues" = mismatch_cnts))

  
  class(result) <- "edit_table"
  return (result)
}