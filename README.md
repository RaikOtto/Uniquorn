# Uniquorn R package

Package to identify cancer cell lines based on their unique somatic mutational fingerprint.

# How to make it work: 

0 Start R, e.g. as R-Studio session

1 Preparation # optional if you have the packages installed and loaded

`install.packages("devtools")`

`library("devtools")`

`source("https://bioconductor.org/biocLite.R")`

`biocLite("VariantAnnotation")`

2 Install the Uniquorn

`install_github("RaikOtto/Uniquorn")`

`library("Uniquorn")`

3 Load example data, here the NCI-60 exome sequencing HT29 Cancer Cell line

`HT29_CL_VCF = paste( system.file("extdata", package="Uniquorn"), "HT29.vcf.gz", sep ="/")`

--> FYI: VCF testfiles for NCI-60 panel e.g. from http://watson.nci.nih.gov/projects/nci60/wes/VCF

4 Run test analysis

`identify_vcf_file( "HT29_cl_vcf"  )`

Under active development!

Contact: raik.otto@hu-berlin.de

Additional information:

You can find suitable DbSNP data for creating your own signatures here http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/
Filename: snp142Common.txt
Place the file in the raw data parser folder
