#' A plotting tool to visualize data from an edit_frame object
#' 
#' @param this an edit_table object
#' @param field character naming member of edit_table to plot. A barplot is produced
#'  for members $AllSites, $RepSites, $mirnaTargetSites, and $VEP. A Venn Diagram
#'  is produced for member $Tissues
#' @param event character naming which column from field to count
#' @param mismatch character indicating which type of mismatch to consider, such as "AtoG".
#'  By default, considers all mismatch types
#' @param group character selecting field to group (color code) by.
#' @param n integer limiting number of events to n. Any number of events beyond
#'  n will be grouped into "other" category
#' @param text_size numeric giving the size of axis and legend labels.
#' @param perc_size numeric giving the size of percent labels above bars.
#' @param perc_round numeric giving decimal places that percent labels round to.
#' @param plot logical where TRUE returns a gg plot object. False returns the modified data.frame from
#'  which the gg object is constructed from. False can be used to combine datasets for "facet plotting".
#' @param x_label character providing label for x axis. Use "none" to disable axis label.
#' @param y_label character providing label for y axis. Use "none" to disable axis label.
#' @param y_range numeric vector providing min and max of y axis range
#' @param denom_all logical. If true, percentage labels will always be
#'  out of the total number of edits found. If false, only the number of edits
#'  in field will be considered.
#' @param legend logical. If false, any legend that would otherwise appear is turned off.
#' @import ggplot2
#' @export
plot.edit_table <- function(this,
                            field = "AllSites",
                            event = "Mismatch",
                            mismatch = "all",
                            group = "Tissue",
                            n = NULL,
                            text_size = 20,
                            perc_size = 10,
                            perc_round = 0,
                            plot = TRUE,
                            x_label = "Mismatch type",
                            y_label = "Number of events",
                            y_range = "auto",
                            denom_all = FALSE,
                            legend = TRUE) {
  
  member <- this[[field]]
  

  if (field == "Tissues")
    stop("Please use tissue_plot() to visulize tissue mismatches with a venn diagram")
  
    
  else {
    
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
    
    # Total num of mismatches, or if mismatch is specified, total number of
    #   mismatches of the specified type.
    if (mismatch != "all")
      tot_cnt <- nrow(member[member$Mismatch == mismatch, ])
    else 
      tot_cnt <- nrow(member)
    
    # If denom_all is true, then tot_cnt should be the total number of
    #   mismatches rather than just those in member
    if (denom_all) {
      if (mismatch != "all")
      tot_cnt <- nrow(this$AllSites[this$AllSites$Mismatch == mismatch, ])
    else 
      tot_cnt <- nrow(this$AllSites)
    }
    
    # Obtain total proportions of each edit across each tissue,
    #   these are used as labels within the plot
    total_prop <- 
      sapply(event_names,
             function(x) {
              cnt <- event_dat[event_dat$Event == x, "Freq"]
              cnt_perc <- round (sum(cnt / tot_cnt) * 100, perc_round) %>%
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
  
    # Lastly, to enable grouping, rename column with name
    #   <group> to "group"
    if (group != "none") {
      i <- which(colnames(event_dat) == group)
      colnames(event_dat)[i] <- "group"
    }
      

    
    if (group != "none") {
      g <- ggplot(event_dat,
                  aes(x = reorder(Event, -Freq),
                      y = Freq,
                      fill = reorder(group, Freq)))
      g <- g + geom_bar(stat = 'identity')
      
    } else {
      g <- ggplot(event_dat,
                  aes(x = reorder(Event, -Freq),
                      y = Freq))
      g <- g + geom_bar(stat = "identity",
                        fill = I("#18453B"),
                        color = I("#18453B"))
    }
    
    g <- g + geom_text(aes(label = Total_prop),
                       na.rm = TRUE,
                       position = "stack",
                       hjust = 0.5,
                       vjust= -0.2,
                       size = perc_size)
    g <- g + ylab(y_label)
    g <- g + xlab(x_label)
    
    if (y_range[1] != "auto")
      g <- g + coord_cartesian(ylim = y_range)
    
    g <- g + theme(axis.title.x = element_text(size = text_size),
                   axis.title.y = element_text(size = text_size),
                   axis.text.x = element_text(angle = 45,
                                              hjust = 1,
                                              size = text_size),
                   axis.text.y = element_text(size = text_size),
                   legend.text = element_text(size = text_size),
                   legend.title = element_text(size = text_size),
                   legend.key.height = unit(3, "line"))
    g <- g + guides(fill = guide_legend(title = "Tissues"))
    
    if (!legend) 
      g <- g + theme(legend.position = "none")
    
    if (y_label == "none")
      g <- g + theme(axis.title.y = element_blank())
    
    if (x_label == "none")
      g <- g + theme(axis.title.x = element_blank())
  }
  if (plot)
    return (g)
  else
    return (cbind(event_dat, field))
}