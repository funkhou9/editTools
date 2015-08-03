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
void rm_search(CharacterMatrix x, std::string rm_file) {
  
  
  std::vector< std::vector<std::string> > mat;
  std::vector< std::vector<std::string> >::iterator vvit;
  std::vector< std::string >::iterator vit;

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
  
  for (vvit = mat.begin(); vvit != mat.end(); vvit++) {
      for (int i = 0; i < x.nrow(); i++) {
        CharacterVector v = x(i, _);
        
        std::vector< std::string > vec(v.size());
        
        for (int i = 0; i < v.size(); i++)
          vec[i] = std::string(v[i]);
        
        if (("chr" + vec[0]) == vvit->at(4) &&
            (std::stol(vec[1]) >= std::stol(vvit->at(5)) && std::stol(vec[1]) <= std::stol(vvit->at(6)))) {
          
          std::string start = vvit->at(5);
          std::string end = vvit->at(6);
          std::string repeat = vvit->at(9);
          std::string family = vvit->at(10);
          
          vec.push_back(start);
          vec.push_back(end);
          vec.push_back(repeat);
          vec.push_back(family);
          
          Rcout << vec << std::endl;
        }
      }
  }
}
  
//   for (int i = 0; i < x.nrow(); i++) {
//     CharacterVector v = x(i, _);
//     
//     std::vector< std::string > vec(v.size());
//     
//     for (int i = 0; i < v.size(); i++)
//       vec[i] = std::string(v[i]);
//     
//     for (vvit = mat.begin(); vvit != mat.end(); vvit++) {
//       
//       if (("chr" + vec[0]) == vvit->at(4) &&
//             (std::stol(vec[1]) >= std::stol(vvit->at(5)) && std::stol(vec[1]) <= std::stol(vvit->at(6)))) {
//         
//         std::string start = vvit->at(5);
//         std::string end = vvit->at(6);
//         std::string repeat = vvit->at(9);
//         std::string family = vvit->at(10);
//         
//         vec.push_back(start);
//         vec.push_back(end);
//         vec.push_back(repeat);
//         vec.push_back(family);
//         
//         Rcout << vec << std::endl;
//       }
//     }
//   }
// }

