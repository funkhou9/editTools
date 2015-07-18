# Filters by genotype
#
# Filter snps that have a particular genotype (homozygous ref, heterozygous, or homozygous var)
# in a particular sample
#
# @param obj vcf object
# @param n sample name or number
# @param geno specify what genotype to filer - 0/0 1/1 or 0/1
# @return filtered vcf object
# @export
gt_filter <- function(obj, n, geno, remove = FALSE) {
  
  # Get genotype
  GT <- samples(obj, n)[[1]]$GT
  
  # Remove or keep selected geno?
  if(remove){
    i <- !(as.character(GT) %in% geno)
  } else i <- as.character(GT) %in% geno
  
  return(obj[i,])
}


# Filters SNPs where samples mismatch
#
# At the moment only works with two sample vcf files
#
# @param obj vcf object
# @return filtered vcf object
# @export
gt_diff_filter <- function(obj) {
  
  # Get 'genotype' call for both DNA and RNA
  mat <- sample_merge(obj, 2, 1)
  
  # Only return those variants 
  i <- as.character(mat[, 1]) != as.character(mat[, 2])
  return(obj[i, ])
}


# Filters out mismatches caused by indels
#
# Doesn't consider more than two alleles?
#
# @param obj vcf object
# @return filtered vcf object
# @export
indel_filter <- function(obj) {
  
  # Both REF and ALT columns need to be just 1 char in length
  bases <- snps(obj)[, c("REF", "ALT")]
  i <- nchar(as.character(bases[, 1])) == 1
  j <- nchar(as.character(bases[, 2])) == 1
  
  return(obj[i & j, ])
}


# Filters on quality
#
# @param obj vcf object
# @param qual quality threshold
# @return filtered vcf object
# @export
qual_filter <- function(obj, qual) {
  
  # A simple quality filter that can be placed on the QUAL column
  i <- snps(obj)[, 6] > qual
  return(obj[i, ])
}


# Filters on RNA coverage - specific to RNA editing work
#
# @param obj vcf object
# @param depth the required depth of RNA in support of a mismatch
# @return filtered vcf object
# @export
edit_depth_filter <- function(obj, depth) {
  
  # Get vector of genotypes
  gt <- as.character(samples(obj, "DNA")[[1]][, "GT"])
  
  # Get RNA total depth and variant depth - ensure numeric
  depedit <- samples(obj, "RNA")[[1]][, c("DP", "DV")]
  depedit <- sapply(depedit, function(x) as.numeric(as.character(x)))
  
  # Convert total depth to reference depth
  depedit <- cbind(DV = depedit[, "DV"],
                   DR = depedit[, "DP"] - depedit[, "DV"])
  
  # building matrix j to index depedit with
  j <- index_geno(gt)
  j <- as.matrix(cbind(seq_along(j), j))
  rownames(j) <- NULL
  colnames(j) <- NULL
  
  # index RNA calls based on DNA
  sampledepth <- depedit[j]
  i <- sampledepth >= depth
  return(obj[i, ])
}


# @param obj vcf object
# @return filtered vcf object
# @export
dp_filter <- function(obj, n, depth, depthv = FALSE) {
  #Minimum depth filter
  #At moment, must specify a single sample and depth filtering
  #   is performed on that sample.
  #By default, overall depth is the filtering criteria. If depthv=T, reads supporting variant become
  #   filtering criteria
  if(depthv) {
    sampledepth <- as.numeric(as.character(samples(obj, n)[[1]][, "DV"]))
  } else sampledepth <- as.numeric(as.character(samples(obj, n)[[1]][, "DP"]))
  i <- sampledepth >= depth
  return(obj[i, ])
}


# @param obj vcf object
# @return filtered vcf object
# @export
hom_filter <- function(obj, n, perc = 95) {
  #Filters snps - a snp must be perc or 1-perc percent homozygous variant based on read counts
  #Must provide a single sample with <n> (order or name)
  depths <- samples(obj, n)[[1]][, c("DP", "DV")]
  pervar <- as.numeric(as.character(depths[, "DV"])) /
              as.numeric(as.character(depths[,"DP"])) * 100
  i <- pervar >= perc | pervar <= (100 - perc)
  return(obj[i, ])
}



