# Younikorn R package

Package to identify cancer cell lines based on their unique somatic mutational fingerprint.

# How to make it work:

1 install.packages("devtools")

2 install_github("RaikOtto/Younikorn")

3 library("Younikorn")

4 download.file( url="http://watson.nci.nih.gov/projects/nci60/wes/VCF/COLO-205.vcf", destfile = "COLO-205.vcf" ) # Downloading a test file 

--> more vcf testfiles from e.g. http://watson.nci.nih.gov/projects/nci60/wes/VCF

5 identify_vcf_file( "COLO-205.vcf" )

Under active development!

Contact: raik.otto@hu-berlin.de
