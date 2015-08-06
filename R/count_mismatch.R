# Counts events for a single tissue sample
#
# @param this_field. Not an edit_table object but a field of an edit_table object, such as
#   this$AllSites or all$RepSites
# @param wname logical. If true, tissue names will be printed on each row.
#   Useful for plotting.
# @param event character naming column to count
# @param trim integer giving the maximum number of events to count. Extras are stored in an
#   "other" category.
# @return data.frame for summarizing results, plotting, etc.
count_mismatch <- function(this_field,
                           wname = FALSE,
                           event = "Mismatch",
                           trim = NULL,
                           mismatch = "all") {

  # Acquire edited tissues
  tiss <- table(this_field[, "Tissue"]) %>%
            sort() %>%
              names()
  

  
  # Construct mismatch table for each
  mismatch_table <-
    lapply(tiss,
           function(x) {
             
             # Two versions of tab. 1st without a mismatch restriction.
             #  The second with a mismatch restriction
             if (mismatch == "all")
               tab <- 
                 this_field[ this_field[, "Tissue"] == x, event ] %>%
                  table() %>%
                    as.data.frame(stringsAsFactors = FALSE)
             else
               tab <- 
                 this_field[ this_field[, "Tissue"] == x &
                                    this_field[, "Mismatch"] == mismatch, event ] %>%
                    table() %>%
                      as.data.frame(stringsAsFactors = FALSE)
               
             
             if (wname) tab <- cbind(tab, x, stringsAsFactors = FALSE) 
             colnames(tab) <- c(event, "Freq", "Tissue")
             
             # Reorder events with most common on top
             tab <- tab [order(tab$Freq, decreasing = TRUE), ]
             
             # Include the proportion of each mismatch
             Prop <- tab$Freq / sum(tab$Freq)
             tab <- cbind(tab, Prop)
             return (tab)
           })
  
  names(mismatch_table) <- tiss


  # First check if trim arg is supplied.
  #   Next check if trim is actually needed
  if (!is.null(trim)) {
    num_events <- sapply(mismatch_table, nrow)
    if (any(trim < num_events)) {
      
      # Find top n events that are common among each sample
      element_list <- 
        lapply(mismatch_table,
               function(x) {
                 x[, event]
               })
      
      common_n_elements <- Reduce(intersect, element_list)[1:trim]
      
      # Shrink mismatch_table to only include common_n_elements
      mismatch_table <-
        lapply(mismatch_table,
               function(x) {
                 # Separate 'common' data from 'other' data
                 tiss_dat <- x
                 tiss_elem <- tiss_dat[, event]
                 idx <- tiss_elem %in% common_n_elements
                 common_dat <- tiss_dat[idx, ]
                 other_dat <- tiss_dat[!idx, ]
                 
                 # Create new 'other' data to paste with 'common' data
                 other_freq <- sum(other_dat$Freq)
                 rbind(common_dat,
                       c("other",
                         other_freq,
                         tiss_dat$Tissue[1],
                         other_freq / (other_freq + common_dat$Freq)))
               })
    }
  }
  return (mismatch_table)
}