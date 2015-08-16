# Read in VEP output intended to be used in edit_table object
#
# @param file character naming VEP output file
read_vep <- function(file) {
  
  # Read all - then prepare to hack away
  tab <- read.table(file,
                    header = TRUE,
                    comment.char = '',
                    sep ='\t',
                    stringsAsFactors = FALSE)
  
  # Split first column into ID, mismatch and tissue type
  tmi <- strsplit(tab[, 1], ':') %>%
            do.call(rbind, .) %>%
              as.data.frame(stringsAsFactors = FALSE)
  colnames(tmi) <- c("Tissue", "Mismatch", "ID")
  
  # Split second column into Chr and Pos - requires two splits
  chr <- strsplit(tab[, 2], ':') %>%
            do.call(rbind, .) %>% 
              as.data.frame(stringsAsFactors = FALSE)
  
  pos <- strsplit(chr[, 2], '-') %>% 
            do.call(rbind, .) %>% 
              as.data.frame(stringsAsFactors = FALSE)
  Pos <- pos[, 1]
  
  cp <- cbind("Chr" = chr[, 1], Pos)
  
  # Get other columns of interest and piece together
  conseq <- tab[, "Consequence"]
  biotype <- tab[, "BIOTYPE"]
  strand <- tab[, "STRAND"]
  gene <- tab[, "SYMBOL"]
  
  # Modify strand to adapt previous '-' and '+' format
  ###
  
  # Modify biotype '-' category to 'unannotated'
  biotype[biotype == '-'] <- "unannotated"
  
  cbind("ID" = tmi[, 3],
        "Chr" = cp[, 1],
        "Pos" = cp[, 2],
        "Strand" = strand,
        "Mismatch" = tmi[, 2],
        "Gene" = gene,
        "Biotype" = biotype,
        "Consequence" = conseq,
        "Tissue" = tmi[, 1]) %>%
    as.data.frame(stringsAsFactors = FALSE)
  
}