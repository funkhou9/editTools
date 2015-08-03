
#' @export
edit_match <- function(editObj, searchObj, stranded=F) {
  # Combs through editObj (output from rnaEditSummary)
  # singleSearch searches searchObj for a single edit
  # returns searchObj row when all "hits" are true
  singleSearch <- function(edit, searchObj, stranded.=stranded) {
    # Looks through searchObj for a single edit site
    # If edit site not in searchObj, there will not
    # be a single position in which hit1, 2 and 3 are all T
    # If this is the case, returns NULL
    
    ### Create variables
    # 'mismatch' vars (atomic values)
    mismatch_pos <- as.numeric(as.character(edit[2]))
    dna <- as.character(unlist(edit[3]))
    rna <- as.character(unlist(edit[4]))
    mismatch <- paste(dna, rna, sep='/')
    mismatch_strand <- as.character(edit[13]) 
    mismatch_chr <- paste("chr", edit[1], sep='')
    
    # target or 'search' vars (vectors)
    search_chr <- as.character(searchObj$chr)
    search_pos1 <- as.numeric(as.character(searchObj$pos1))
    search_pos2 <- as.numeric(as.character(searchObj$pos2))
    search_strand <- as.character(searchObj$strand)
    
    # Find 'hits' in target
    hit1 <- mismatch_chr == search_chr
    hit2 <- mismatch_pos >= search_pos1
    hit3 <- mismatch_pos <= search_pos2
    hitmat <- data.frame(hit1, hit2, hit3)
    
    if (stranded.) {
      hit4 <- mismatch_strand == search_strand
      hitmat <- cbind(hitmat, hit4)
      if (max(rowSums(hitmat)) < 4) {
        return(NULL)
      } else {
        result <- cbind(mismatch_pos, mismatch, searchObj[hit1 & hit2 & hit3 & hit4,])  
        return (result)
      }
    } else if (max(rowSums(hitmat)) < 3) {
      return(NULL)
    } else {
      result <- cbind(mismatch_pos, mismatch, searchObj[hit1 & hit2 & hit3,])
      return(result)
    }
  }
  
  # Storing just edit calls (not summary)
  allEdits <- editObj$Calls
  
  # apply singleSearch across all edits
  resultList <- apply(allEdits, 1, singleSearch, searchObj)
  
  # stores result in data.frame
  result <- do.call(rbind.data.frame, resultList)
  return(result)
}