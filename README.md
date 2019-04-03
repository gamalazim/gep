# gep
**gep** is a package built in C with R interface. Legacy fortran fspack libraries are included. The package can be used to solve massive systems of sparse mixed model equations. *gep* offers control over memory and takes extra precautions to save solutions on disk - not only in R workspace.    

**To install *gep* directly from github:**

```r
devtools::install_github("gamalazim/gep")
library(gep)
```