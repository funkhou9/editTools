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
  
  member <- this[[field]]
  
  # Names of each tissue in this object. Sorted from lowest edits to highest
  tiss_names <- table(member$Tissue) %>%
    sort() %>%
      names()
  
  # Each type of mismatch present in the dataset
  mismatch_names <- unique(member[, "Mismatch" ]) 
  
  # Total num of mimsatch
  tot_cnt <- nrow(member)
  
  # Mismatch types present in each tissue
  types <- sapply(tiss_names,
                  function(x) {
                    member[member$Tissue == x, "Mismatch"] %>%
                      unique()
                  })
  
  # Number of mismatch types in each tissue
  type_counts <- sapply(types, length)
  

  if (field == "Tissues") {

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
    
    # Obtain each tissue summary
    mismatches <- lapply(tiss_names,
                         count_mismatch,
                         this[[field]],
                         wname = TRUE)
    
    mismatches <- do.call(rbind, mismatches)
    names(mismatches)[3] <- "Tissue"
 
    # Obtain total proportions of each edit across each tissue,
    #   these are used as labels within the plot
    total_prop <- 
      sapply(mismatch_names,
             function(x) {
              cnt <- mismatches[mismatches$Mismatch == x, "Freq"]
              round (sum(cnt / tot_cnt) * 100, 2) %>%
                paste('%', sep = '')
           })
    
    # Obtain index for each mismatch type - where is it represented last
    #   in the mismatches object? (Necessary for adding total_prop to mismatches)
    type_idx <- 
      sapply(mismatch_names,
            function(x) {
              which(mismatches$Mismatch == x) %>%
                max()
           })
    
    # Add total_prop to each row of mismatches, appropriately
    mismatches[type_idx, "Total_prop"] <- total_prop
    
    
    ggplot(mismatches,
           aes(x = reorder(Mismatch, -Freq),
               y = Freq,
               fill = Tissue)) +
      geom_bar(stat = 'identity') +
      geom_text(aes(label = Total_prop),
                na.rm = TRUE,
                position = "stack",
                hjust = 0.5,
                vjust= -0.2) +
      ylab("Number of events") +
      xlab("Type of mismatch") +
      theme(axis.title.x = element_text(size=15),
            axis.title.y = element_text(size=15),
            axis.text.x = element_text(angle = 45,
                                       hjust = 1,
                                       size=15),
            axis.text.y = element_text(size=15),
            strip.text.x = element_text(size=15),
            legend.text = element_text(size=15),
            legend.title = element_text(size=15)) +
      guides(fill = guide_legend(title = "Tissues"))
  }
}