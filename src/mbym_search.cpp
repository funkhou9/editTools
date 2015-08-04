#include <Rcpp.h>
using namespace Rcpp;

#include <fstream>


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
void mbym_search(CharacterMatrix& x, std::string rm_file) {
  
  
  std::vector< std::vector<std::string> > mat;
  std::vector< std::vector<std::string> >::iterator vvit_inner;
  std::vector< std::vector<std::string> >::iterator vvit_outer;

  std::vector< std::vector<std::string> > query_mat;
  
  long l_num = 0;
  
  std::string line;
  std::ifstream rm(rm_file);
  
  while (getline(rm, line)) {
    
    l_num ++;
    
    if (l_num >= 4) {
      std::vector< std::string > rm_line = parse_r(line);
      mat.push_back(rm_line);
    }
  }
  
  for (int i = 0; i < x.nrow(); i++) {
    CharacterVector v = x(i, _);
    std::vector< std::string > vect_q(v.size());
    
    for (int i = 0; i < v.size(); i++) {
      vect_q[i] = std::string(v[i]);
    }
    
    query_mat.push_back(vect_q);
  }

  for (vvit_inner = query_mat.begin(); vvit_inner != query_mat.end(); vvit_inner++) {
    for (vvit_outer = mat.begin(); vvit_outer != mat.end(); vvit_outer++) {
    
      if (("chr" + vvit_inner->at(0)) == vvit_outer->at(4) &&
          (std::stol(vvit_inner->at(1)) >= std::stol(vvit_outer->at(5)) && std::stol(vvit_inner->at(1)) <= std::stol(vvit_outer->at(6)))) {
        
        std::string start = vvit_outer->at(5);
        std::string end = vvit_outer->at(6);
        std::string repeat = vvit_outer->at(9);
        std::string family = vvit_outer->at(10);
        
        (*vvit_inner).push_back(start);
        (*vvit_inner).push_back(end);
        (*vvit_inner).push_back(repeat);
        (*vvit_inner).push_back(family);
        
        Rcout << (*vvit_inner) << std::endl;
      }
    }
  }
}

