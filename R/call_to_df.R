# Used in read_vcf to allocate n numbers of samples
#
# @param x the columns of a vcf object which contain sample information
# @return parsed sample information
call_to_df <- function(x, names=NULL) {
  
  parse <- strsplit(as.character(x), ':')
  parsed_allele <- data.frame(matrix(unlist(parse),
                                     nrow = length(parse),
                                     byrow = TRUE),
                              stringsAsFactors = FALSE)
  names(parsed_allele) <- names
  return(parsed_allele)
}