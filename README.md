# Younikorn R package

Package to identify cancer cell lines based on their unique somatic mutational fingerprint.

# How to make it work:

1 Preparation # optional if you have the packages installed and loaded

`install.packages("devtools")`

library("devtools")

source("https://bioconductor.org/biocLite.R")

biocLite("VariantAnnotation")

3 install_github("RaikOtto/Younikorn")

4 library("Younikorn")

5 download.file( url="http://watson.nci.nih.gov/projects/nci60/wes/VCF/COLO-205.vcf", destfile = "COLO-205.vcf" ) # Downloading a test file 

--> more vcf testfiles from e.g. http://watson.nci.nih.gov/projects/nci60/wes/VCF

6 identify_vcf_file( "COLO-205.vcf" )

Under active development!

Contact: raik.otto@hu-berlin.de
