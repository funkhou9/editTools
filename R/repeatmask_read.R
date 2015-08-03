
#' @export
repeatmask_read <- function(file, skip = 3, chrs = c(1:18, "X")) {
  # Stores repeatmasker data as data frame
  # Takes only certain columns from repeatmasker (.out) data
  # Retains only chromosomes of interest in argument <chrs>
  rmask <- read.table(file, header=F, skip=skip, sep='', nrow=100)
  classes <- sapply(rmask, class)
  rmask <- read.table(file, header=F, skip=skip, sep='', colClasses=classes)
  chromosomes <- paste("chr", chrs, sep='')
  rmask <- rmask[,c(5,6,7,10,11)]
  names(rmask) <- c("chr", "pos1", "pos2", "repeat", "family")
  rmask <- rmask[rmask$chr %in% chromosomes,]
  return(rmask)  
}