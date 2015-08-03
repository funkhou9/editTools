
#' @export
repeatmask_read <- function(file, skip = 3, chrs = c(1:18, "X")) {
  
  # Read sample of file to peek at column classes
  rmask <- read.table(file,
                      header = FALSE,
                      skip = skip,
                      sep = '',
                      nrow = 100)
  
  classes <- sapply(rmask, class)
  
  # Re-read with col classes
  rmask <- read.table(file,
                      header = FALSE,
                      skip = skip,
                      sep='',
                      colClasses = classes)
  
  chromosomes <- paste("chr", chrs, sep='')
  
  # Only grabbing position begin/end, repeat and familiy columns
  rmask <- rmask[, c(5, 6, 7, 10, 11)]
  
  names(rmask) <- c("chr", "pos1", "pos2", "repeat", "family")
  
  # Only observing those chromosomes specified
  rmask <- rmask[rmask$chr %in% chromosomes, ]
  return(rmask)
}