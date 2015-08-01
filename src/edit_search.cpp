#include "Variant.h"

// @export
//' @useDynLib editTools
//' @importFrom Rcpp sourceCpp
// [[Rcpp::export]]
void edit_search(std::string file,
                char strand,
                CharacterVector names,
                bool ex_indel,
                int geno_dp,
                int geno_hom,
                int edit_dp,
                int lh)
{
  
  std::string line;
  std::ifstream vcf1(file);
  std::vector< std::string > header;
  
  std::vector< std::string > names_vec(names.size());
  for (int i = 0; i < names.size(); i++) {  
    names_vec[i] = std::string(names[i]);  
  }
  
  // Additional criteria for filtering
  // Requires homozygous genotypes
  std::vector< std::string > genos;
  genos.push_back("0/0");
  genos.push_back("1/1");
  
  // For each line, check if a ## header line and exclude
  while (getline(vcf1, line)) {
    if (line.empty() || (line.find("##")) == 0)
      continue;
    
    // When encounitering # header line, store field names
    if(line.find('#') == 0) {
      header = parse_v(line);
      continue;
    }
    
    // Initialize Variant object, flagged with strand information
    Variant Var(line, strand);
    
    // Parse the rest of the line and add to list of rna samples
    std::vector< std::string > line_vec = parse_v(line);
  
    // Initialize Rna objects with one of two ways depending on if
    //  names is provided
    if (names_vec.empty()) {
      for (int i = 10; i < line_vec.size(); i++) {
        Rna r(line_vec[i], header[i]);
        Var.add_rna(r);
      }
    } else {
      for (int i = 0; i < names_vec.size(); i++) {
        Rna r(line_vec[i + 10], names_vec[i]);
        Var.add_rna(r);
      }
    }
    
    // Flag Rna objects for evidence for editing
    Var.gt_diff_filter();
    Var.edit_depth_filter(edit_dp);
    Var.likelihood_filter(lh);
    
    // If Variant passes all filters, call genotypes for each sample
    //  and print
    if (Var.indel_filter() &&
        // Var.qual_filter(qual) &&
        Var.gt_filter(genos) &&
        Var.hom_filter(geno_hom) &&
        Var.dp_filter(geno_dp) &&
        Var.contains_edit()) {
      
      Var.call_samples();
      Rcout << Var;
    }
  }
}
