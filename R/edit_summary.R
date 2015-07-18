# Scans vcf class for evidence of RNA editing
# 
# @param obj A vcf object
# @param qual An integer specifiying the minimum variant QUAL
# @param ex.indel logical indicating whether to exclude indels from the scan
# @param geno.dp integer specifying the minimum genotype depth
# @param geno.hom integer ranging from 0 to 1 specifiying the proportion of homozygosity
#  the genotype must exhibit
# @param edit.dp integer specifying the minimum depth required for evidence of
#  RNA editing
# @param summary logical indicating whether a summary of found edits is returned
#  or a subsetted vcf object
# @return Either a subsetted vcf object, only containing probable edited loci
#  or an object of class "edit_summary"
edit_summary <- function(obj,
                         qual_,
                         ex.indel_,
                         geno.dp_,
                         geno.hom_,
                         edit.dp_,
                         summary_) {
  
  obj_edit <-
    gt_filter(obj, "DNA", c("0/0", "1/1")) %>%
      hom_filter("DNA", geno.hom_) %>%
        gt_diff_filter() %>%
          indel_filter() %>%
            qual_filter(qual_) %>%
              dp_filter("DNA", geno.dp_) %>%
                edit_depth_filter(edit.dp_)
                  
  if (summary_) {
    result <- call_geno(obj_edit, length(obj_edit$Samples))
    return(result)
    
  } else {
    return(obj_edit)
  }
}