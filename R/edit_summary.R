#' Scans vcf class for evidence of RNA editing
#' 
#' @param obj A vcf object
#' @param qual An integer specifiying the minimum variant QUAL
#' @param ex.indel logical indicating whether to exclude indels from the scan
#' @param geno.dp integer specifying the minimum genotype depth
#' @param geno.hom integer ranging from 0 to 1 specifiying the proportion of homozygosity
#'  the genotype must exhibit
#' @param edit.dp integer specifying the minimum depth required for evidence of
#'  RNA editing
#' @param summary logical indicating whether a summary of found edits is returned
#'  or a subsetted vcf object
#' @return Either a subsetted vcf object, only containing probable edited loci
#'  or an object of class "edit_summary"
#' @import magrittr
#' @export
edit_summary <- function(obj,
                         qual = 10,
                         ex.indel = TRUE,
                         geno.dp = 10,
                         geno.hom = 95,
                         edit.dp = 5,
                         summary = TRUE) {
  
  obj_edit <-
    gt_filter(obj, "DNA", c("0/0", "1/1")) %>%
      hom_filter("DNA", geno.hom) %>%
        gt_diff_filter() %>%
          indel_filter() %>%
            qual_filter(qual) %>%
              dp_filter("DNA", geno.dp) %>%
                edit_depth_filter(edit.dp)
                  
  if (summary) {
    result <- call_geno(obj_edit, length(obj_edit$Samples))
    class(result) <- "edit_summary"
    return(result)
  } else {
    return(obj_edit)
  }
}