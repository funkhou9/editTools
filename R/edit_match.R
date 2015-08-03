
#' @export
edit_match <- function(editObj, searchObj, stranded=F) {
  
  

  
  # Storing just edit calls (not summary)
  allEdits <- editObj$Calls
  
  # apply singleSearch across all edits
  resultList <- apply(allEdits, 1, single_match, searchObj)
  
  # stores result in data.frame
  result <- do.call(rbind.data.frame, resultList)
  return(result)
}