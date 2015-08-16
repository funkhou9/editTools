#' Adds Variant Effect Predictor information to an edit_table object
#' 
#' @param this an edit_table object
#' @param vep_file character giving path to txt file with VEP results
#' @return a new edit_table object with added $VEP field
#' @export
add_vep <- function(this, vep_file) {
  
  vep <- read_vep(vep_file)
  
  new_result <- append(this,
                       list("VEP" = vep),
                       1)
  class(new_result) <- "edit_table"
  return (new_result)
}