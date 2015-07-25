#include <Rcpp.h>
using namespace Rcpp;

#include <fstream>

#include "Variant.h"

//' @export
//' @useDynLib editTools
//' @importFrom Rcpp sourceCpp
// [[Rcpp::export]]
void edit_search(std::string file,
                char strand,
                long qual = 10,
                bool ex_indel = true,
                int geno_dp = 10,
                int geno_hom = 95,
                int edit_dp = 5)
{
  
  std::string line;
  std::ifstream vcf1(file);
  
  // Additional criteria for filtering
  // Requires homozygous genotypes
  std::vector<std::string> genos;
  genos.push_back("0/0");
  genos.push_back("1/1");
  
  // For each line, check if a header line
  while (getline(vcf1, line)) {
    if (line.empty() || (line.find('#')) == 0)
      continue;
    
    // Initialize Variant object bool flagged as "plus strand" variants
    Variant Var(line, strand);
    
    if (Var.indel_filter() &&
        Var.qual_filter(qual) &&
        Var.gt_filter(genos) &&
        Var.hom_filter(geno_hom) &&
        Var.dp_filter(geno_dp) &&
        Var.gt_diff_filter() &&
        Var.edit_depth_filter(edit_dp)) {
      
      Rcout << Var;
    }
  }
}
