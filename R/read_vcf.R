#' Read in .vcf files
#' 
#' @param filename The filename of the .vcf file
#' @param NS The number of samples included in the .vcf file
#' @param names Vector of names that will be used to reference each sample. 
#' Specify in the order that the samples appear in the file.
#' @return An object of class vcf
#' @export
read_vcf <- function(filename, NS, names=NULL) {
  bcf <- read.table(filename, header=F, sep='\t')
  snps <- bcf[, 1:8]
  colnames(snps) <- c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO")
  format <- unlist(strsplit(as.character(bcf[1,9]), ":"))
  samples <- apply(bcf[, 10:ncol(bcf)], 2, function(x) call.to.df(x, names=format))
  names(samples) <- names
  result <- list(SNPs = snps, Samples = samples)
  class(result) <- "vcf"
  return(result)
}