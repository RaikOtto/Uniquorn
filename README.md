# Uniquorn R package

Package to identify cancer cell lines based on their unique somatic mutational fingerprint.

# How to make it work: 

0 Start R, e.g. as R-Studio session

1 Preparation # optional if you have the packages installed and loaded

`install.packages("devtools")`

`library("devtools")`

`source("https://bioconductor.org/biocLite.R")`

`biocLite("VariantAnnotation")`

2 Installing and loading the Uniquorn

`install_github("RaikOtto/Uniquorn")`

`library("Uniquorn")`

3 Run test analysis

`data("HT29_cl_vcf",package="Uniquorn")`

--> FYI: VCF testfiles for NCI-60 panel e.g. from http://watson.nci.nih.gov/projects/nci60/wes/VCF

`identify_vcf_file( "HT29_cl_vcf"  )`

Under active development!

Contact: raik.otto@hu-berlin.de
