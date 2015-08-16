# editTools - Detection and visualization of RNA editing data from VCF files

### Installation
editTools contains compiled code and relies on the Rcpp package  
If using GNU version 4.7 or later, set c++11 with
```R
Sys.setenv("PKG_CXXFLAGS" = "-std=c++11")
```
Otherwise
```R
Sys.setenv("PKG_CXXFLAGS" = "-std=c++0x")
```

Install with
```R
devtools::install_github("funkhou9/editTools")
```

See ?find_edits for primary usage


