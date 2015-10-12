/********************************************************************** 
 * RNA editing detector for VCF data
 *
 * Goals:
 *  - Ecapulate the complexity of 'variants' presented in VCF files
 *  - Provide a class to be used by editTools::edit_search()
 *    a C++ function to be used in the R interpreter
 *  
 *  
 * Author: Scott Funkhouser <funkhou9@msu.edu>
 **********************************************************************/

#include <Rcpp.h>
using namespace Rcpp;

#include <fstream>

/**********************************************************
 * Global parse_v() - Used to split vcf lines into vectors
 *  Overloaded (2 flavors)
 *  1. Provide a string
 *  2. Provide a string and a delimiter
 *  
 *  See definitions at bottom of file
 **********************************************************/

// Intended for general space separated fields
std::vector< std::string > parse_v(const std::string& line);


// Intended for special colon (or other) separated fields
std::vector< std::string > parse_v(const std::string& line,
                                   char sep);

char delim_samp = ':';
char delim_pl = ',';
char delim_info = ';';
char delim_equals = '=';



class Rna
{
  
  /**************************************************
   * Representing a single RNA locus in a VCF file
   * 
   * rna_gt - coded genotype for a variant
   * rna_pl - genotype likelihood
   * rna_dp - total sequencing depth for a variant
   * rna_dv - sequencing depth in support of variant
   * call - base call
   * diff_flag - true if call differs from corresponding
   *  DNA sample
   * depth_flag - if edit depth meets criteria
   * likelihood_flag - if samples genotype likelihoods
   *  meet criteria
   * edit_frac - proportion of reads that support edit
   **************************************************/

public:
  std::string tissue_name;
  std::string rna_gt;
  std::vector< std::string > rna_pl;
  double rna_dp;
  double rna_dv;
  double edit_dp;
  std::string call;
  bool diff_flag;
  bool depth_flag;
  // bool likelihood_flag;
  double edit_frac;
  double sb;
  bool sb_flag;

  
public:
  Rna(const std::string& line, std::string tissue_name)
  {
    this->diff_flag = false;
    this->depth_flag = false;
    // this->likelihood_flag = false;
    
    std::vector< std::string > rna_call = parse_v(line, delim_samp);
    
    this->rna_gt = rna_call[0];
    this->rna_pl = parse_v(rna_call[1], delim_pl);
    this->rna_dp = std::stod(rna_call[2]);
    this->rna_dv = std::stod(rna_call[3]);
    this->sb = std::stod(rna_call[4]);
    this->sb_flag = false;
    
    this->tissue_name = tissue_name;
  }
};



class Variant
{
  
  /************************************************************
   * VCF file produces these fields when options 
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
   ************************************************************/
  
  std::string chrom;
  unsigned long pos;
  std::string ref;
  std::string alt;
  long qual;
  std::string dna_gt;
  std::vector< std::string > dna_pl;
  double dna_dp;
  double dna_dv;
  std::list< Rna > rna_list;
  std::string call;
  bool geno_likelihood_flag;
  int ave_mq;
  
  // Additional variable for strand ID, 
  //  either '+' or '-'.
  char strand;
  
public:
   
  // Initialize with a line from a vcf file and strand ID
  Variant(const std::string& line, char& strand)
  {
    this->strand = strand;
    this->geno_likelihood_flag = false;

    // Parse whole line into 'general' fields
    std::vector< std::string > gen_set = parse_v(line);
    
    // Distribute gen_set elements
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
    this->dna_pl = parse_v(dna_call[1], delim_pl);
    this->dna_dp = std::stod(dna_call[2]);
    this->dna_dv = std::stod(dna_call[3]);
    
    // Search info field for helpful tags
    std::vector< std::string > info = parse_v(gen_set[7],
                                              delim_info);
    
    this->ave_mq = std::stoi(parse_v(info.back(), delim_equals).at(1));
  }
  
  
  // Add an RNA sample
  void add_rna(Rna& r)
  {
    this->rna_list.push_back(r);
  }
  

  
  /* **************************************************
   * Filters for RNA editing detection
   * 
   * Each method flags true if Variant object
   *  provides a piece of evidence for RNA editing
   *****************************************************/
  
  
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
  
  // Detects if Variant possesses an Rna object in rna_list where its genotype doesn't
  //  match the genomic sample
  void gt_diff_filter()
  {
    for (std::list<Rna>::iterator it = rna_list.begin(); it != rna_list.end(); it++) {
      if (dna_gt != it->rna_gt)
        it->diff_flag = true;
    }
  }
    
  // Detects if Variant possesses an Rna object in rna_list where the depth of sequence
  //  supporting RNA editing is at least the depth specified
  void edit_depth_filter(int& depth)
  {
    
    if (dna_gt == "0/0")
      for (std::list<Rna>::iterator it = rna_list.begin(); it != rna_list.end(); it++) {
        it->edit_dp = it->rna_dv;
        it->edit_frac = it->edit_dp / it->rna_dp;
        if (it->edit_dp >= depth)
          it->depth_flag = true;
      }
        
    else
      for (std::list<Rna>::iterator it = rna_list.begin(); it != rna_list.end(); it++) {
        it->edit_dp = it->rna_dp - it->rna_dv;
        it->edit_frac = it->edit_dp / it->rna_dp;
        if (it->edit_dp >= depth)
          it->depth_flag = true;
      }
  }

  // Sufficient likelihood of RNA call? Second most likely genotype call must have a phred
  //  scaled genotype likelihood greater than l.
//   void likelihood_filter(int& l)
//   {
//     
//     // Inspect likelihoods for genome
//     int count = 0;
//     
//     for (std::vector< std::string >::const_iterator itp = dna_pl.begin(); itp != dna_pl.end(); itp++) {
//       int likeli = std::stoi(*itp);
//       if (likeli >= l)
//         count ++;
//     }
//     
//     if (count == 2)
//       geno_likelihood_flag = true;
//     
//     // Inspect likelihoods for each rna sample
//     for (std::list<Rna>::iterator it = rna_list.begin(); it != rna_list.end(); it++) {
//       int count_rna = 0;
//       
//       for (std::vector< std::string >::const_iterator itp = it->rna_pl.begin(); itp != it->rna_pl.end(); itp++) {
//         int likeli = std::stoi(*itp);
//         if (likeli >= l)
//           count_rna ++;
//       }
//       
//       if (count_rna == 2)
//         it->likelihood_flag = true;
//     }
//   }
  
  void sb_flag(int bias)
  {
    for (std::list<Rna>::iterator it = rna_list.begin(); it != rna_list.end(); it++) {
      if (it->sb <= bias)
        it->sb_flag = true;
    }
  }
  
  // If the Variant object has at least one Rna object in its list that meets both
  //  criteria for editing
  bool contains_edit()
  {
    for (std::list<Rna>::iterator it = rna_list.begin(); it != rna_list.end(); it++) {
      if (it->depth_flag && it->diff_flag && it->sb_flag)
        return true;
    }
    return false;
  }
  
  /*********************************************************** 
   *  Methods to call genomic and RNA samples
   ***********************************************************/
   
  
  // Changes base calls to their reverse compliment
  void reverse_strand()
  {
    if (ref == "A") ref.assign("T");
    else if (ref == "T") ref.assign("A");
    else if (ref == "C") ref.assign("G");
    else if (ref == "G") ref.assign("C");

    if (alt == "A") alt.assign("T");
    else if (alt == "T") alt.assign("A");
    else if (alt == "C") alt.assign("G");
    else if (alt == "G") alt.assign("C");
  }
  

  // Provides RNA call by returning either
  //  REF or ALT, depending on rna_gt and strand.
  // If rna_gt is heterozygous, reports base call
  //  that is supportive of RNA editing.
  void call_samples()
  {
    
    if (strand == '-')
      reverse_strand();
      
    // Call DNA sample
    if (dna_gt == "0/0") 
      this->call = ref;
    else
      this->call = alt;
    
    // Call each RNA sample
    for (std::list<Rna>::iterator it = rna_list.begin(); it != rna_list.end(); it++) {
      if (it->rna_gt == "0/1") {
        
        if (dna_gt == "0/0") {
          it->call = alt;
        } else {
          it->call = ref;
        }
    
      } else if (it->rna_gt == "0/0") {
        it->call = ref;
        
      } else {
        it->call = alt;
      }
    }
  }
  

  // Variant objects send certain members and method results to stdout
  //  Can be used with cout or Rcout. Output intended to be redirected.
  friend std::ostream& operator<<(std::ostream& os, Variant& var)
  {
    for (std::list<Rna>::iterator it = var.rna_list.begin(); it != var.rna_list.end(); it++) {
      if (it->depth_flag && it->diff_flag && it->sb_flag)
        os << var.chrom << '\t' << var.pos << '\t' << var.strand <<
          '\t' <<  var.call << "to" << it->call << '\t' << var.dna_dp << '\t' << var.dna_dv << '\t' <<
            it->rna_dp << '\t' << it->edit_dp << '\t' << it->edit_frac << '\t' << it->sb << '\t' <<
              var.ave_mq << '\t' << it->tissue_name << std::endl;
    }
    
    return os;
  }
};



/************************************************************ 
 * Global definitions
 ************************************************************/


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
