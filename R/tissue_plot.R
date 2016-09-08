#' Plotting tool to produce venn diagrams showing off tissue specific mismatch counts
#'
#' @param this an edit_table object
#' @param field character naming edit_table field (either $AllSites,
#'  $RepSites, or $mirnaTargetSites) from which to obtain counts from
#' @param fill_colors character vector specifying colors for each tissue
#' @param cex_labels numeric specifying size of tissue labels
#' @param cex_counts numeric specifying size of mismatch counts
#' @param mismatch character. "all" counts all mismatches. Specific mismatches
#'  can be specified, for example "AtoG".
#' @return NULL
#' @export
tissue_plot <- function(this,
                        field = "AllSites",
                        fill_colors = "green",
                        cex_labels = 2,
                        cex_counts = 3.5,
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
  
  venn <- VennDiagram::venn.diagram(positions,
                                    filename = NULL,
                                    fill = fill_colors,
                                    alpha = 0.7,
                                    cat.cex = cex_labels,
                                    cex = cex_counts,
                                    cat.fontfamily = "Helvetica",
                                    sub.fontfamily = "Helvetica")

  plot.new()
  grid::grid.draw(venn)
  
}
