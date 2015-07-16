#' Getter for SNP field
#
#' @param obj a vcf object
#' @return SNPs field of a vcf object
#' @export
snps <- function(obj) obj$SNPs

#' Getter for Samples field
#
#' @param obj a vcf object
#' @param n either a numeric or character identifying which
#'  sample to extract
#' @return Samples field (a list of samples) or if n provided,
#'  a data.frame with information for a single sample.
#' @export
samples <- function(obj, n)  obj$Samples[n]

#' Subsetting method for vcf class
#' 
#' @param obj a vcf obect
#' @param p numeric to filter SNPs
#' @param n numeric to filter Samples
#' @return a subsetted vcf object
#' @export
`[.vcf` <- function(obj, p, n) {
  #Method for `[` that will be used in other functions to filter and subset
  L <- list(SNPs = snps(obj)[p, ], 
            Samples = lapply(samples(obj, n), function(x){x[p, ]}))
  class(L) <- "vcf"
  return(L)
}

