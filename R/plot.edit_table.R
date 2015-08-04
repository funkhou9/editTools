#' A plotting tool to visualize DNA RNA mismatch counts from an edit_frame
#' 
#' @import ggplot2
#' @export
plot.edit_table <- function(this) {
  
  # Rearrange levels of both Mismatch and Tissue in order to present
  #   from highest to lowest
  ggplot(this,
         aes(reorder(Mismatch, Mismatch, function(x) -length(x)))) +
    geom_bar(aes(fill = reorder(Tissue, Tissue, function(x) - length(x)))) +
    ylab("Number of events") +
    xlab("Type of mismatch") +
    theme(axis.title.x = element_text(size=15),
          axis.title.y = element_text(size=15),
          axis.text.x = element_text(angle = 45, hjust = 1, size=15),
          axis.text.y = element_text(size=15),
          strip.text.x = element_text(size=15),
          legend.text = element_text(size=15),
          legend.title = element_text(size=15))
  
  
}