#' Plotting tool to produce venn diagrams showing off tissue specific mismatch counts
#'
#' @param this an edit_table object
#' @param field character naming edit_table field (either $AllSites,
#'  $RepSites, or $mirnaTargetSites) from which to obtain counts from
#' @param fill_colors character vector specifying colors for each tissue
#' @return NULL
#' @import VennDiagram
#' @import grid
#' @export
tissue_plot <- function(this,
                        field = "AllSites",
                        fill_colors = "green",
                        mismatch = "all") {
  
  member <- this[[field]]
  
  tiss_names <- member$Tissue %>%
                  table() %>%
                    sort() %>% 
                      names()
  
  
  positions <- 
    lapply(tiss_names,
           function(x) {
             if (mismatch == "all") {
               paste(member[member$Tissue == x, "Chr"],
                     member[member$Tissue == x, "Pos"],
                     member[member$Tissue == x, "Strand"],
                     member[member$Tissue == x, "Mismatch"])
             } else {
               paste(member[member$Tissue == x & member$Mismatch == mismatch, "Chr"],
                     member[member$Tissue == x & member$Mismatch == mismatch, "Pos"],
                     member[member$Tissue == x & member$Mismatch == mismatch, "Strand"],
                     member[member$Tissue == x & member$Mismatch == mismatch, "Mismatch"])
             }
           })
  
  names(positions) <- tiss_names
  
  venn <- venn.diagram(positions,
                       filename = NULL,
                       fill = fill_colors,
                       alpha = 0.7,
                       cat.cex = 2,
                       cex = 3.5,
                       cat.fontfamily = "Helvetica",
                       sub.fontfamily = "Helvetica")

  plot.new()
  grid.draw(venn)
  
}
