# editTools
### Detection and visualization of RNA editing data from VCF files

### Installation
editTools contains compiled code and relies on the Rcpp package and c++11.    
If using GNU version 4.7 or later, specify c++11 with:
```R
Sys.setenv("PKG_CXXFLAGS" = "-std=c++11")
```
Otherwise, use:
```R
Sys.setenv("PKG_CXXFLAGS" = "-std=c++0x")
```

Then install the editTools package with:
```R
devtools::install_github("funkhou9/editTools")
```

See ?find_edits for primary usage


