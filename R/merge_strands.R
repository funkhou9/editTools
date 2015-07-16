#' Merges edit_summary objects from different strands
#' 
#' @param plus_obj an edit_summary object - generated from RNA sequencing reads
#'  that mapped to the PLUS strand
#' @param minus_obj an edit_summary object - generaged from RNA sequencing reads
#'  that mapped ot the MINUS strand
#' @return a new edit_summary object that is the result of properly merging
#'  both plus_obj and minus_obj
#' @export
merge_strands <- function(plus_obj, minus_obj) {

  # First correct minus_obj
  minus <- correct_sequence(minus_obj)
  
  # Obtain Calls field from each
  minus_calls <- minus$Calls
  plus_calls <- plus_obj$Calls
  
  # Produce vectors of '-' or '+' used for a new column
  #   in output
  min_strand <- rep('-', nrow(minus_calls))
  plus_strand <- rep('+', nrow(plus_calls))
  
  # Add in strand column
  minus_calls$Strand <- min_strand
  plus_calls$Strand <- plus_strand
  
  # plus and minus results stacked on top of each other - reorder
  callsresult <- rbind(plus_calls,
                       minus_calls)
  callsresult <- callsresult[order(callsresult[, 1], callsresult[, 2]), ]
  
  # Modifying corresponding Summary field
  minus_summary <- minus$Summary
  plus_summary <- plus_obj$Summary
  rownames(minus_summary) <- paste(minus_summary[,"DNA"],
                                   minus_summary[,"RNA"],
                                   sep = 'to')
  
  rownames(plus_summary) <- paste(plus_summary[,"DNA"],
                                  plus_summary[,"RNA"],
                                  sep = 'to')
  
  # Reorder summaries so that they have the same order
  minus_summary <- minus_summary[rownames(plus_summary), ]
  
  
  Freq <- minus_summary$Freq + plus_summary$Freq
  prob <- Freq / sum(Freq)
  summresult <- cbind(minus_summary[ ,1:2], Freq, prob)
  
  result <- list("Calls" = callsresult,
                 "Summary" = summresult)
  class(result) <- "edit_summary"
  return (result)
}
