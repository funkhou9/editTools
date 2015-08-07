#' Add miRNA target data to an edit_table object
#'
#' @param this a edit_table object
#' @param mirna_file character string naming formatted miranda ouptut file
#' @return a newly formatted edit_table object with miRNA target information
#' @export
add_mirna <- function(this, mirna_file) {
  
  mod_result <-
    capture.output(mbym_search(as.matrix(this$AllSites), mirna_file, 2, 3, 4, 3, 4, 0, 1, TRUE, 5)) %>%
      sapply(function(x) strsplit(x, split = '\t')) %>%
        do.call(rbind, .) %>%
          as.data.frame(stringsAsFactors = FALSE)
  
  # Format output - replace sequence depth data with repeatmasker data
  mod_result <- mod_result[, c(1, 2, 3, 4, 5, 9, 10, 11, 12, 8)]
  rownames(mod_result) <- NULL  
  colnames(mod_result) <- c("ID",
                            "Chr",
                            "Pos",
                            "Strand",
                            "Mismatch",
                            "Target_start",
                            "Target_end",
                            "miRNA",
                            "Target_gene",
                            "Tissue")
  
  new_result <- append(this,
                       list("mirnaTargetSites" = mod_result),
                       1)
  
  class(new_result) <- "edit_table"
  
  return (new_result)
}