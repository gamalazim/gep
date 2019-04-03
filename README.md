# gep
**gep** is a package built in C with R interface. Legacy fortran fspak libraries are included. The package can be used to solve massive systems of sparse mixed model equations. *gep* offers control over memory and takes extra precautions to save solutions on disk - not only in R workspace.    

**To install *gep* directly from github:**

```r
devtools::install_github("gamalazim/gep")
library(gep)
```
Expect a few compiler warnings due to the legacy fortran code in fspak libraries! Feel free to adjust the code yourself, otherwise they are 'benign' warnings! 
