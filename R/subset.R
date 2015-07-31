#' Subset editing candidates from an edit_table object
#' 
#' Can subset based on such criteria as tissue type or DNA call
#' 
#' If subsetted both by DNA and RNA call (mismatch type) then a
#'  data.frame is returned rather than an edit_table object
#' 
#' @param this an edit_table object
#' @param dna character indicating desired DNA call
#' @param rna character indicating desired RNA call
#' @param tissue character indicating desired tissue sample
#' @export
#' @return a subsetted edit_table
subset.edit_table <- function(this,
                              dna = NULL,
                              rna = NULL,
                              tissue = NULL) {
  this$Edits[, "DNA"] == dna
  this$Edits[, "RNA"] == rna
  this$Edits[, "Tissue"] == tissue
}