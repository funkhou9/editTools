# Produces $counts field - counts of each DNA/RNA mismatch type
#
# Produces $counts for a single sample. Vectorized in find_edits()
count_mismatch <- function(n, edits) {

  # Construct mismatch table for tissue n
  mismatch_table <-
    edits [edits[, "Tissue"] == n, "Mismatch"] %>%
      table() %>%
        data.frame()
    
  colnames(mismatch_table) <- c("Mismatch", "Freq")
  
  # Reorder so most common mismatches are displayed on top
  mismatch_table <- mismatch_table [order(mismatch_table$Freq, decreasing = TRUE), ]
  
  # Include the proportion of each mismatch
  mismatch_table$Freq <- as.character(mismatch_table$Freq) %>% as.numeric()
  Prop <- mismatch_table$Freq / sum(mismatch_table$Freq)
  
  mismatch_table <- cbind(mismatch_table, Prop)
  
  return (mismatch_table)
}