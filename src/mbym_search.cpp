#include <Rcpp.h>
using namespace Rcpp;

#include <fstream>
#include <string>
#include <vector>


std::vector< std::string > parse_r(const std::string& line)
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

std::ostream& operator<<(std::ostream& os, std::vector< std::string >& field)
{
  for (std::vector< std::string >::const_iterator it = field.begin(); it != field.end(); it++) {
    os << *it << '\t';
  }
  return os;
}


// [[Rcpp::export]]
void mbym_search(CharacterMatrix& x,
                 std::string rm_file,
                 int s_chr,
                 int s_start,
                 int s_end,
                 int item_1,
                 int item_2,
                 int item_3,
                 int item_4,
                 bool stranded = false,
                 int s_strand = 5)
{
  
  long l_num = 0;
  std::vector< std::vector<std::string> > mat, query_mat;
  std::vector< std::vector<std::string> >::iterator vvit_query, imin, imax, imid;
  std::string line;
  std::ifstream rm(rm_file);
  std::string pos1;
  std::string pos2;
  std::string repeat;
  std::string family;
  
  // Read file into 2D vector of strings
  while (getline(rm, line)) {
    l_num ++;
    
    if (l_num >= 4) {
      std::vector< std::string > rm_line = parse_r(line);
      mat.push_back(rm_line);
    }
  }
  
  // Convert R char matrix to 2D vector of strings
  for (int i = 0; i < x.nrow(); i++) {
    CharacterVector v = x(i, _);
    std::vector< std::string > vect_q(v.size());
    
    for (int i = 0; i < v.size(); i++) {
      vect_q[i] = std::string(v[i]);
    }
    query_mat.push_back(vect_q);
  }
  
  // multi-tiered binary search-like approach
  // First search for chromosome, then position for each row of query matrix
  for (vvit_query = query_mat.begin(); vvit_query != query_mat.end(); vvit_query++) {
    // Reset subject iterators
    imin = mat.begin();
    imax = mat.end();
    while (std::distance(imin, imax) != 1) {
      int mid = std::distance(imin, imax) / 2;
      imid = imin;
      std::advance(imid, mid);
      
      // If chromosome and position match...
      if ("chr" + vvit_query->at(1) == imid->at(s_chr) &&
         (std::stol(vvit_query->at(2)) >= std::stol(imid->at(s_start))) &&
         (std::stol(vvit_query->at(2)) <= std::stol(imid->at(s_end)))) {
        // Grab info from subject,
        (*vvit_query).push_back(imid->at(item_1));
        (*vvit_query).push_back(imid->at(item_2));
        (*vvit_query).push_back(imid->at(item_3));
        (*vvit_query).push_back(imid->at(item_4));
        
        if (stranded) {
          if (vvit_query->at(3) == imid->at(s_strand)) {
            Rcout << (*vvit_query) << std::endl;
            break;
          }
        } else {
          // Then print and break
          Rcout << (*vvit_query) << std::endl;
          break;
        }
        
      } else if ("chr" + vvit_query->at(1) > imid->at(s_chr) ||
                (std::stol(vvit_query->at(2)) > std::stol(imid->at(s_end)))) {
        imin = imid++;
        
      } else {
        imax = imid--;
      }
    }
  }
}


