/* RNA editing detector for VCF data
 * 
 * Goals:
 *  - Ecapulate the complexity of 'variants' presented in VCF files
 *  - Provide a class to be used by editTools::edit_search()
 *    a C++ function to be used in the R interpreter
 *  
 * Note: Preprocessor directives in edit_search.cpp for easier use with
 *  Rcpp
 *  
 * Author: Scott Funkhouser <funkhou9@msu.edu>
 */

class Variant
{
  
  /* VCF file produces these fields when options
   *  used are -O v -m -v with two samples (DNA and RNA in this case) 
   * 
   * chrom = chromosome where variant was found
   * pos = chromosomal position where variant was found
   * id = some sort of identifier... ?
   * ref = base or bases present at this position in the reference genome
   * alt = alternative base(s) discovered at this position
   * qual = assessment of confidence in variant call. Higher is better.
   * filter = processed by other software... ?
   * info = ';' delimited sequence of additional information
   * format = ';' delimited sequence listing how dna_call and rna_call
   *  should be read 
   */
  
  std::string chrom;
  unsigned long pos;
  std::string ref;
  std::string alt;
  long qual;
  std::string dna_gt;
  std::string dna_pl;
  double dna_dp;
  double dna_dv;
  std::string rna_gt;
  std::string rna_pl;
  double rna_dp;
  double rna_dv;

  // Additional variable for strand ID, 
  //  either '+' or '-'.
  char strand;
  
 public:
  
  // Initialize with a line from a vcf file and strand ID
  Variant(const std::string& line, char& strand)
  {
    
    this->strand = strand;
    
    // Prepare vcf line for tokenizing
    std::istringstream iss(line);
    std::vector< std::string > attr_set;
    
    // Deposit tokens into 'general' mem. variables
    while (iss) {
      std::string attr;
      iss >> attr;
      
      attr_set.push_back(attr);
    }
    
    // Distribute attr_set elements
    this->chrom = attr_set[0];
    this->pos = stol(attr_set[1]);
    this->ref = attr_set[3];
    this->alt = attr_set[4];
    this->qual = stol(attr_set[5]);

    // Prepare DNA sample for tokenizing
    std::istringstream dna_st(attr_set[9]);
    std::string dna_tag;
    std::vector<std::string> dna_call;
    
    // Deposit tokens into dna mem. variables
    while (getline(dna_st, dna_tag, ':')) {
      dna_call.push_back(dna_tag);
    }
    
    // Distribute dna_call elements
    this->dna_gt = dna_call[0];
    this->dna_pl = dna_call[1];
    this->dna_dp = stod(dna_call[2]);
    this->dna_dv = stod(dna_call[3]);
    
    // Prepare RNA sample for tokenizing
    std::istringstream rna_st(attr_set[10]);
    std::string rna_tag;
    std::vector<std::string> rna_call;
    
    // Deposit tokens into rna mem. variables
    while (getline(rna_st, rna_tag, ':')) {
      rna_call.push_back(rna_tag);
    }
    
    // Distribute rna_call elements
    this->rna_gt = rna_call[0];
    this->rna_pl = rna_call[1];
    this->rna_dp = stod(rna_call[2]);
    this->rna_dv = stod(rna_call[3]);
  }
  
  
  
  
  /* Filters for RNA editing detection
   * 
   * Each method flags true if Variant object
   *  provides a piece of evidence for RNA editing
   */
  
  
  // Detects if Variant has a single char in both
  //  REF and ALT fields
  bool indel_filter()
  {
    return (ref.length() == alt.length() == 1);
  }
  
  // Detects if Variant meets a quality threshold
  bool qual_filter(long& qual_thresh)
  {
    return (qual >= qual_thresh);
  }
  
  // Detects if Variant has specified genotype(s)
  bool gt_filter(std::vector< std::string >& geno)
  {
    for (int i = 0; i < geno.size(); i++) {
      if (geno[i] == dna_gt)
        return true;
    }
    return false;
  }
  
  // Detects if Variant is homozygous according
  //  to a specified percentage of sequencing reads
  bool hom_filter(int& perc)
  {
    double perc_var = (dna_dv / dna_dp) * 100;
    return (perc_var >= perc || perc_var <= (100 - perc));
  }
  
  // Detects if Variant has a DNA sample has sufficient depth
  bool dp_filter(int& depth)
  {
    return (dna_dp >= depth);
  }
  
  // Detects if Variant has differing DNA and RNA calls
  bool gt_diff_filter()
  {
    return (dna_gt != rna_gt);
  }
  
  // Detects if Variant has a sufficient number of RNA reads
  //  that are in support of editing
  bool edit_depth_filter(int& depth)
  {
    if (dna_gt == "0/0")
      return(rna_dv >= depth);
    else
      return(rna_dp - rna_dv >= depth);
  }
  
  
  
  /* 
   *  
   *  Methods to choose between and change REF
   *    and ALT fields
   *  
   *
   */
   
  
  // Changes base calls to their reverse compliment
  void reverse_strand()
  {
    if (ref == "A") ref = "T";
    if (ref == "T") ref = "A";
    if (ref == "C") ref = "G";
    if (ref == "G") ref = "C";
    
    if (alt == "A") alt = "T";
    if (alt == "T") alt = "A";
    if (alt == "C") alt = "G";
    if (alt == "G") alt = "C";
  }
  
  // Provides DNA call by returning either
  //  REF or ALT, depending on dna_gt and strand
  std::string call_dna() {
    
    if (strand == '-') {
      reverse_strand();
    }
    
    if (dna_gt == "0/0")
      return ref;
    else
      return alt;
  }
   
  // Provides RNA call by returning either
  //  REF or ALT, depending on rna_gt and strand.
  // If rna_gt is heterozygous, reports base call
  //  that is supportive of RNA editing.
  std::string call_rna()
  {
    if (strand == '-') {
      reverse_strand();
    }
    
    if (rna_gt == "0/1") {
      
      if (dna_gt == "0/0") {
        return alt;
      } else {
        return ref;
      }
      
    } else if (rna_gt == "0/0") {
      return ref;
      
    } else {
      return alt;
    }
  }
  
  // Variant objects send certain members  and method results to stdout
  //  Can be used with cout or Rcout. Output intended to be redirected.
  friend std::ostream& operator<<(std::ostream& os, Variant& var)
  {
    os << var.chrom << '\t' << var.pos << '\t' << var.call_dna() << '\t' << var.call_rna() <<
      '\t' << var.dna_dp << '\t' << var.dna_dv << '\t' <<
        var.rna_dp << '\t' << var.rna_dv << std::endl;
    
    return os;
  }
};

