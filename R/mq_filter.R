#' Filters potential RNA editing sites by applying filtering on average mapping quality (MQ)
#' 
#' Calculates P(X >= x | N, p(e)) where x is the number of reads which provide support for
#' editing for sample i at site j, N is the total depth for sample i at site j,
#' and p(e) is the probability of error for site j, derived from the average mapping quality.
#'
#' @param this an edit_table object 
#' @export
mq_filter <- function(this) {
  
  # Obtained average quality for each site and convert to probability
  ave_mq <- as.numeric(this$AllSites$Ave_MQ)
  pe <- 10^-(ave_mq / 10)
  
  # Total rna depths and mismatch depths
  rna_dp <- as.numeric(this$AllSites$RNA_depth)
  edit_dp <- as.numeric(this$AllSites$RNA_mismatch_depth)
  
  # Estimate conditional probability P(X >= x | N, p(e)) ~ Binom(N, p(e))
  1 - pbinom(edit_dp, rna_dp, pe)
}