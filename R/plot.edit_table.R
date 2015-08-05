#' A plotting tool to visualize DNA RNA mismatch counts from an edit_frame
#' 
#' @param this an edit_table object
#' @param character naming member of edit_table to plot. A barplot is produced
#'  for members $AllSites, $RepSites, $mirnaTargetSites, and $VEP. A Venn Diagram
#'  is produced for member $Tissues
#' @import ggplot2
#' @import grid
#' @import VennDiagram
#' @export
plot.edit_table <- function(this, field = "AllSites") {
  
  if (field == "Tissues") {
    tiss_names <- table(this$AllSites$Tissue) %>%
                    sort(decreasing = TRUE) %>%
                      names()
    
    positions <- 
      lapply(tiss_names,
             function(x) {
               paste(this$AllSites[this$AllSites$Tissue == x, "Chr"],
                     this$AllSites[this$AllSites$Tissue == x, "Pos"],
                     this$AllSites[this$AllSites$Tissue == x, "Strand"],
                     this$AllSites[this$AllSites$Tissue == x, "Mismatch"])
             })
    
    names(positions) <- tiss_names
    
    venn <- venn.diagram(positions,
                        filename = NULL,
                        fill = rainbow(length(tiss_names)),
                        alpha = 0.65,
                        cat.cex = 2,
                        cex = 2.5,
                        cat.fontfamily = "Arial")
    
    plot.new()
    grid.draw(venn)
    
  } else {
    
    # Rearrange levels of both Mismatch and Tissue in order to present
    #   from highest to lowest
    ggplot(this[[field]],
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
            legend.title = element_text(size=15)) +
      guides(fill = guide_legend(title = "Tissues"))
    
  }
  

  
  
}