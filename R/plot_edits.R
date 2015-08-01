#' A plotting tool to visualize DNA RNA mismatch counts from an edit_frame
#' 
#' @import ggplot2
#' @export
plot_edits <- function(this) {
  
  ggplot(this, aes(Mismatch), fill = Tissue) +
    geom_bar()
  
}