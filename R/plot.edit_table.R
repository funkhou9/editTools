#' A plotting tool to visualize data from an edit_frame object
#' 
#' @param this an edit_table object
#' @param field character naming member of edit_table to plot. A barplot is produced
#'  for members $AllSites, $RepSites, $mirnaTargetSites, and $VEP. A Venn Diagram
#'  is produced for member $Tissues
#' @param event character naming which column from field to count
#' @param mismatch character indicating which type of mismatch to consider, such as "AtoG".
#'  By default, considers all mismatch types
#' @param n integer limiting number of events to n. Any number of events beyond
#'  n will be grouped into "other" category
#' @import ggplot2
#' @import grid
#' @import VennDiagram
#' @export
plot.edit_table <- function(this,
                            field = "AllSites",
                            event = "Mismatch",
                            mismatch = "all",
                            n = NULL) {
  
  member <- this[[field]]
  

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
    
    # Obtain each tissue summary with added args to count_match()
    event_dat <- count_mismatch(member,
                                wname = TRUE,
                                event = event,
                                trim = n,
                                mismatch = mismatch)
    
    
    
    # Join list and ensure proper formatting needed for plotting
    event_dat <- do.call(rbind, event_dat)
    event_dat[, c("Freq")] <- as.numeric(event_dat[, c("Freq")])
    event_dat[, c("Prop")] <- as.numeric(event_dat[, c("Prop")])
    colnames(event_dat)[1] <- "Event"
  
    
    # Names of each tissue in this object. Sorted from lowest edits to highest
    tiss_names <- unique(event_dat$Tissue)
    
    # Each type of mismatch present in the dataset
    event_names <- unique(event_dat$Event) 
    
    # Total num of events in all of data field
    tot_cnt <- nrow(member)
    

    # Obtain total proportions of each edit across each tissue,
    #   these are used as labels within the plot
    total_prop <- 
      sapply(event_names,
             function(x) {
              cnt <- event_dat[event_dat$Event == x, "Freq"]
              cnt_perc <- round (sum(cnt / tot_cnt) * 100, 2) %>%
                            paste('%', sep = '')
              return (cnt_perc)
             })
    
    
    # Obtain index for each mismatch type - where is it represented last
    #   in the event_dat object? (Necessary for adding total_prop to event_dat)
    type_idx <- 
      sapply(event_names,
             function(x) {
              which(event_dat$Event == x) %>%
                max()
            })
    
    # Add total_prop to each row of event_dat, appropriately
    event_dat[type_idx, "Total_prop"] <- total_prop
  

    g_plot <-
      ggplot(event_dat,
             aes(x = reorder(Event, -Freq),
                 y = Freq,
                 fill = reorder(Tissue, Freq))) +
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
  
  return (g_plot)
}