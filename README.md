# Younikorn R package

Package to identify cancer cell lines based on their unique somatic mutational fingerprint.

# How to make it work: 

0 Start R, e.g. as R-Studio session

1 Preparation # optional if you have the packages installed and loaded

`install.packages("devtools")`

`library("devtools")`

`source("https://bioconductor.org/biocLite.R")`

`biocLite("VariantAnnotation")`

2 Procuring test data

`download.file( url="http://watson.nci.nih.gov/projects/nci60/wes/VCF/COLO-205.vcf", destfile = "COLO-205.vcf" )` # Downloading a test file

--> more vcf testfiles from e.g. http://watson.nci.nih.gov/projects/nci60/wes/VCF

3 Installing and loading the Younikorn

`install_github("RaikOtto/Younikorn")`

`library("Younikorn")`

4 Run test analysis

`identify_vcf_file( "COLO-205.vcf" )`

Under active development!

Contact: raik.otto@hu-berlin.de
