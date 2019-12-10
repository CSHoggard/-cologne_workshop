### Momocs

if(!require("Momocs")) install.packages('Momocs', repos='http://cran.us.r-project.org') 

### tidyverse

if(!require("tidyverse")) install.packages('tidyverse', repos='http://cran.us.r-project.org')

### cowplot
 
if(!require("cowplot")) install.packages('cowplot', repos='http://cran.us.r-project.org')

### geomorph (BETA edition)*

install.packages("devtools")

devtools::install_github("geomorphR/geomorph", ref = "Develop", build_vignettes = TRUE)

### *Due to issues with rgl (necessary for digitisation and visualisation), please use the developer copy and not the CRAN version.

### ggtree (via BiocManager)

install.packages("BiocManager")
library(BiocManager)
install("ggtree")