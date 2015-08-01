# Produces $counts field - counts of each DNA/RNA mismatch type
#
# Produces $counts for a single sample. Vectorized in find_edits()
count_mismatch <- function(n, edit_frame) {

  # Construct mismatch table for tissue n
  mismatch_table <-
    edit_frame [edit_frame$Tissue == n, c("DNA", "RNA")] %>%
      table() %>%
        data.frame()
    
  # Remove 0 counts (those A/A, C/C, G/G, T/T non-mismatches)
  mismatch_table <- mismatch_table [mismatch_table$Freq != 0, ]
  
  # Replace DNA and RNA columns with 'mismatch' column
  mismatch_results <- 
    paste(mismatch_table$DNA, mismatch_table$RNA, sep = 'to') %>%
      cbind(mismatch_table$Freq) %>%
        as.data.frame()
  
  colnames(mismatch_results) <- c("Mismatch", "Freq")
  
  # Reorder so most common mismatches are displayed on top
  mismatch_results <- mismatch_results [order(mismatch_results$Freq, decreasing = TRUE), ]
  
  # Include the proportion of each mismatch
  mismatch_results$Freq <- as.character(mismatch_results$Freq) %>% as.numeric()
  Prop <- mismatch_results$Freq / sum(mismatch_results$Freq)
  
  mismatch_results <- cbind(mismatch_results, Prop)
  
  return (mismatch_results)
}