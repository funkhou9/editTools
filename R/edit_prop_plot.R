#' Plots distributions of mismatch fractions
#' 
#' Mismatch fractions are defined as the proportion of
#'  reads that map to a locus that support a mismatch out of all
#'  the reads that map to that locus
#'  
#' @param this an edit_table object
#' @param field character specifiying edit_table field. Options are
#'  "AllSites" to consider all mismatches, "Repsites" to consider
#'  only mismatches within repetitive regions, and "mirnaTargetSites"
#'  to consider only mismatches in putative 3'UTR mirna Target sites
#' @param use.nonAtoG logical. Do you want to group all non AtoG mismatches together?
#' @param line_size numeric specifying size of lines in plot
#' @return a ggplot functional
#' @export
edit_prop_plot <- function(this,
                           field = "AllSites",
                           use.nonAtoG = TRUE,
                           text_size = 20,
                           line_size = 2) {
  
  member <- this[[field]]
  
  # Get freq_dat from AllSites so that it can be used in any field
  freq_dat <-
    this$AllSites[, c("ID", "Mismatch", "RNA_edit_frac", "Tissue")]
  
  # Reassemble new df that contains sites from member
  new_df <- 
    freq_dat[freq_dat$ID %in% member$ID, ]
  
  # Ensure edit_frac is numeric
  new_df$RNA_edit_frac <- as.numeric(new_df$RNA_edit_frac)
  
  # Check logicals
  if (use.nonAtoG) 
    new_df[new_df$Mismatch != "AtoG", "Mismatch"] <- "non-AtoG"
  
  g <- ggplot(new_df,
              aes(RNA_edit_frac, color = Mismatch))
  
  g <- g + geom_freqpoly(fill = "transparent",
                         size = line_size)
  g <- g + ylab("Count")
  g <- g + xlab("Mismatch Proportion")
  g <- g + theme(axis.title.x = element_text(size = text_size),
                 axis.title.y = element_text(size = text_size),
                 axis.text.y = element_text(size = text_size),
                 axis.text.x = element_text(size = text_size),
                 strip.text.x = element_text(size = text_size),
                 legend.text = element_text(size = text_size),
                 legend.title = element_text(size = text_size))
  return (g)
}