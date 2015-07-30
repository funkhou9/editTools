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
  
  # Reorder so most common mismatches are displayed on top
  mismatch_table <- mismatch_table [order(mismatch_table$Freq, decreasing = TRUE), ]
  
  return (mismatch_table)
}