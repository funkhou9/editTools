#' Add repeatmasker hits to an edit_table object
#'
#' @param this a edit_table object
#' @param rm_file character string naming repeatmasker file
#' @return a newly formatted edit_table object with repeatmasker information
#' @export
add_repeatmask <- function(this, rm_file) {
  
  mod_result <-
    capture.output(mbym_search(as.matrix(this$AllSites), rm_file, 4, 5, 6, 5, 6, 9, 10)) %>%
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
                            "Repeat_start",
                            "Repeat_end",
                            "Element",
                            "Family",
                            "Tissue")
  
  new_result <- append(this,
                      list("RepSites" = mod_result),
                      1)
  
  class(new_result) <- "edit_table"
  
  return (new_result)
}
