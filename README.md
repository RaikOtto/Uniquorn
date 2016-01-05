# Uniquorn R package

Package to identify cancer cell lines based on their weighted mutational fingerprint.

# How to make it work: 

Start R

# 1 Preparation 

 optional if you have the packages installed and loaded

`install.packages("devtools")`

`library("devtools")`

`source("https://bioconductor.org/biocLite.R")`

# 2 Install the Uniquorn

`install_github("RaikOtto/Uniquorn")`

`library("Uniquorn")`

# 3 Test run

 Here the NCI-60 exome sequenced HT29 Cancer Cell line

`HT29_CL_VCF = paste( system.file("extdata", package="Uniquorn"), "HT29.vcf.gz", sep ="/")`

--> FYI: VCF testfiles for NCI-60 panel e.g. from http://watson.nci.nih.gov/projects/nci60/wes/VCF

# 4 Run test analysis

`identify_vcf_file( HT29_CL_VCF  )`
You will find an output file with the name of your input VCF + the ending _uniquorn_identification.tab next to the input VCF file if you did not specify the output file path.

# 5 Add CCLE and CoSMIC CLP CL data

Please add the CCLE and CoSMIC CLP Cancer Cell Line (CL) data manually due to legal regulations! Else only the vanilla 62 CellMiner CLs will be available for identification. You can however manually add custom CLs.

Download into e.g. '~/Downloads' the files 

'CosmicCLP_MutantExport.tsv'

from http://cancer.sanger.ac.uk/cell_lines/download

'CCLE_hybrid_capture1650_hg19_allVariants_2012.05.07.maf'

http://www.broadinstitute.org/ccle/data/browseData

R

`setwd('~/Downloads')`

`initiate_canonical_databases()`

Contact: raik.otto@hu-berlin.de


