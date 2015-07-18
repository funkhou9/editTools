# Merging tool to bring together genotypes (homo var, het, or homo ref) from
# different samples together.
#
# @param object a vcf object
# @param number of samples in vcf file
# @param col which columns would you like to merge from each sample? examples:1,2,3 or "GT", "PL"
# @return a merged sample object
sample_merge <- function(obj, NS, col) {
  
  # Grab specified columns
  cols <- lapply(samples(obj, 1:NS),
                 function(x) { data.frame(x[, col]) })
  
  # Merge
  call <- do.call("cbind", cols)
  return(call)
}