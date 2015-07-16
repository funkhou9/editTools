# Corrects a edit_summary object from a 'minus' strand vcf file
# 
# VCF files always report variant locations according to what is on the TOP
#  strand. Results from minus strand transcripts need to be corrected
# 
# @param editObj
# @param strand
# @return list resembling vcf object with corrected base calls
correct_sequence <- function(editObj, strand = "minus") {
  if(strand == "minus") {
    calls <- editObj$Calls
    summ <- editObj$Summary
    
    calls[, "DNA"] <- 
      unname(sapply(as.character(calls[, "DNA"]),
                    switch, 'A'='T','T'='A','G'='C','C'='G'))
    
    calls[, "RNA"] <-
      unname(sapply(as.character(calls[, "RNA"]),
                    switch, 'A'='T','T'='A','G'='C','C'='G'))
    
    summ[, "DNA"] <-
      unname(sapply(as.character(summ[, "DNA"]),
                    switch, 'A'='T','T'='A','G'='C','C'='G'))
    
    summ[, "RNA"] <-
      unname(sapply(as.character(summ[, "RNA"]),
                    switch, 'A'='T','T'='A','G'='C','C'='G'))
    return(list("Calls" = calls,
                "Summary" = summ))
  }
}