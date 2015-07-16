#' Read in .vcf files
#' 
#' Initializes a vcf object for use in editTools methods
#' Currently, only works with 2-sample VCF files
#' 
#' @param filename The filename of the .vcf file
#' @param NS The number of samples included in the .vcf file
#' @param names character vector of names that will be used to reference each sample. 
#' Specify in the order that the samples appear in the file.
#' @return An object of class vcf
#' @export
read_vcf <- function(filename, NS, names = NULL) {
  
  # Only handles smaller files (-v option to only report variants)
  bcf <- read.table(filename, header=F, sep='\t')
  
  # 'SNP info' (not sample specific) stored in first 8 columns
  snps <- bcf[, 1:8]
  colnames(snps) <- c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO")
  
  # format contains information for how to read sample information
  format <- unlist(strsplit(as.character(bcf[1, 9]), ":"))
  
  # separate samples into a list and provide format names for columns
  samples <- apply(bcf[, 10:ncol(bcf)], 2,
                   function(x) call_to_df(x, names = format))
  
  # name each sample according to provided names
  names(samples) <- names
  
  # vcf class - a 2 element list. The first contains a
  #   df of SNP information. The second is a list of 2
  #   (number of samples), each a data.frame with
  #   sample specific information formatted to "format"
  result <- list(SNPs = snps, Samples = samples)
  class(result) <- "vcf"
  return(result)
}