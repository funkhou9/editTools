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



/* Global parse_v() - Used to split vcf lines into vectors
 *  Overloaded (2 flavors)
 *  1. Provide a string
 *  2. Provide a string and a delimiter
 *  
 *  See implementations below
 */

// Intended for general space separated fields
std::vector< std::string > parse_v(const std::string& line);


// Intended for special colon separated fields
std::vector< std::string > parse_v(const std::string& line,
                                   char sep);



/* Global << overloadings 
 * To print a vector of strings
 *  
 *  See implementations below
 */

std::ostream& operator<<(std::ostream& os, std::vector< std::string >& field);


class Rna
{
  /* Representing a single RNA sample in a VCF file
   * 
   * rna_gt coded genotype for a variant
   * rna_pl ...
   * rna_dp total sequencing depth for a variant
   * rna_dv sequencing depth in support of variant
   */

public:     
  std::string rna_gt;
  std::string rna_pl;
  double rna_dp;
  double rna_dv;
  
public:
  Rna(const std::string& line)
  {
    char delim_samp = ':';
    std::vector< std::string > rna_call = parse_v(line, delim_samp);
    
    this->rna_gt = rna_call[0];
    this->rna_pl = rna_call[1];
    this->rna_dp = std::stod(rna_call[2]);
    this->rna_dv = std::stod(rna_call[3]);
  }
};


std::ostream& operator<<(std::ostream& os, std::list< Rna >& field)
{
  for (std::list<Rna>::iterator it = field.begin(); it != field.end(); it++) {
    os << it->rna_dp << '\t' << it->rna_dv << '\t';
  }
  return os;
}


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
  std::list< Rna > rna_list;
  
  // Additional variable for strand ID, 
  //  either '+' or '-'.
  char strand;

public:
  std::vector< std::string > calls;

public:
   
  // Initialize with a line from a vcf file and strand ID
  Variant(const std::string& line, char& strand)
  {
    char delim_samp = ':';
    
    this->strand = strand;
    
    // Parse whole line into 'general' fields
    std::vector< std::string > gen_set = parse_v(line);
    
    // Distribute attr_set elements
    this->chrom = gen_set[0];
    this->pos = std::stol(gen_set[1]);
    this->ref = gen_set[3];
    this->alt = gen_set[4];
    this->qual = std::stol(gen_set[5]);

    // Parse DNA information
    std::vector< std::string > dna_call = parse_v(gen_set[9],
                                                  delim_samp);
    
    // Distribute dna_call elements
    this->dna_gt = dna_call[0];
    this->dna_pl = dna_call[1];
    this->dna_dp = std::stod(dna_call[2]);
    this->dna_dv = std::stod(dna_call[3]);
  }
  
  
  // Add an RNA sample
  void add_rna(Rna& r)
  {
    this->rna_list.insert(rna_list.begin(), r);
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
    for (std::list<Rna>::iterator it = rna_list.begin(); it != rna_list.end(); it++) {
      if (dna_gt != it->rna_gt)
        return true;
    }
    return false;
  }
    
  // Detects if Variant has a sufficient number of RNA reads
  //  that are in support of editing
  bool edit_depth_filter(int& depth)
  {
    if (dna_gt == "0/0")
      for (std::list<Rna>::iterator it = rna_list.begin(); it != rna_list.end(); it++) {
        if (it->rna_dv >= depth)
          return true;
      }
        
    else
      for (std::list<Rna>::iterator it = rna_list.begin(); it != rna_list.end(); it++) {
        if (it->rna_dp - it->rna_dv >= depth)
          return true;
      }

    return false;
  }

  
  /* 
   *  
   *  Methods to choose between and change REF and ALT fields
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
  

  // Provides RNA call by returning either
  //  REF or ALT, depending on rna_gt and strand.
  // If rna_gt is heterozygous, reports base call
  //  that is supportive of RNA editing.
  void call_samples()
  {

    if (strand == '-') {
      reverse_strand();
    }
    
    // Call DNA sample
    if (dna_gt == "0/0") 
      this->calls.push_back(ref);
    else
      this->calls.push_back(alt);
    
    // Call each RNA sample
    for (std::list<Rna>::iterator it = rna_list.begin(); it != rna_list.end(); it++) {
      if (it->rna_gt == "0/1") {
        
        if (dna_gt == "0/0") {
          this->calls.push_back(alt);
        } else {
          this->calls.push_back(ref);
        }
    
      } else if (it->rna_gt == "0/0") {
        this->calls.push_back(ref);
        
      } else {
        this->calls.push_back(alt);
      }
    }
  }
  

  // Variant objects send certain members and method results to stdout
  //  Can be used with cout or Rcout. Output intended to be redirected.
  friend std::ostream& operator<<(std::ostream& os, Variant& var)
  {
    os << var.chrom << '\t' << var.pos << '\t' << var.calls <<
      var.dna_dp << '\t' << var.dna_dv << '\t' << 
        var.rna_list << std::endl;
    
    return os;
  }
};



/* 
 * Global definitions
 */

std::ostream& operator<<(std::ostream& os, std::vector< std::string >& field)
{
  for (std::vector< std::string >::const_iterator it = field.begin(); it != field.end(); it++) {
    os << *it << '\t';
  }
  return os;
}


std::vector< std::string > parse_v(const std::string& line)
{
  std::istringstream iss(line);
  std::string attr;
  std::vector< std::string > attr_set;
  
  // Deposit tokens into vector
  while (iss) {
    iss >> attr;
    attr_set.push_back(attr);
  }
  return attr_set;
}


std::vector< std::string > parse_v(const std::string& line,
                                   char sep)
{
  std::istringstream iss(line);
  std::string attr;
  std::vector< std::string > attr_set;
  
  // Deposit tokens into a vector, remove sep 
  while (getline(iss, attr, sep)) {
    attr_set.push_back(attr);
  }
  return attr_set;
}  
