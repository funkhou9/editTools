


# @export
single_match <- function(edit, search_obj, stranded_) {

  # Mismatch attributes
  mismatch <- as.character(edit[4])
  mismatch_pos <- as.numeric(as.character(edit[2]))
  mismatch_strand <- as.character(edit[3]) 
  mismatch_chr <- paste("chr", edit[1], sep = '')
  
  # Search obj attributes
  search_chr <- as.character(search_obj$chr)
  search_pos1 <- as.numeric(as.character(search_obj$pos1))
  search_pos2 <- as.numeric(as.character(search_obj$pos2))
  search_strand <- as.character(search_obj$strand)
  
  # Find 'hits' in target
  hit1 <- mismatch_chr == search_chr
  hit2 <- mismatch_pos >= search_pos1
  hit3 <- mismatch_pos <= search_pos2
  hitmat <- data.frame(hit1, hit2, hit3)
  
  if (stranded_) {
    hit4 <- mismatch_strand == search_strand
    hitmat <- cbind(hitmat, hit4)
    if (max(rowSums(hitmat)) < 4) {
      return(NULL)
    } else {
      
      # Save results - sub out DP and DV columns with repeat information
      result <- cbind(edit[c(1, 2, 3, 4)],
                      search_obj[hit1 & hit2 & hit3 & hit4, ],
                      edit[9])  
      return (result)
    }
  } else if (max(rowSums(hitmat)) < 3) {
    return(NULL)
  } else {
    result <- cbind(mismatch_pos,
                    mismatch,
                    search_obj[hit1 & hit2 & hit3, ])
    return(result)
  }
}