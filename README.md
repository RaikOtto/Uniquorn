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

Note that some systems require the prior command `options(unzip = 'internal')` to install from github 

# 3 Test run

Here the NCI-60 exome sequenced HT29 Cancer Cell line

`library("Uniquorn")`

`data("HT29_vcf_fingerprint")`

`ident_result = identify_vcf_file( HT29_vcf_file  )`

`head( ident_result )` will show a table with potential identification candidate, how many mutations overall and weighted of the training set have been found and if any training samples have surpased the identification threshold.

The HT29 cancer cell line vcf file which contains the somatic mutations of the HT-29 cancer cell line is taken from http://watson.nci.nih.gov/projects/nci60/wes/VCF
The Watson repository contains the same cancer cell line samples as Uniquorn's default training dataset, the original CellMiner panel, available @ http://discover.nci.nih.gov/cellminer/. 
However, the cancer cell line vcf sample was differentely filtered and different algorithms have been used to predict mutations and variations. Therefore, it was used as example set to show the difficulties associated with an identification of cancer cell line sample: Some mutations are found only in the query version (here HT-29 from Watson) and some only in the training dataset (here HT-29 from CellMiner original panel). Moreover, some mutations of the query cancer cell line map to different cancer cell lines in the CellMiner training database what can lead to an incorrect identification of a cancer cell line. 

Therefore, a robust yet sensitive cancer cell line identification algorithm is required.

You will find a file with the ending '_uniquorn_identification.tab' next to the input VCF file if you did not specify the output file path.

# 4 Add CCLE and CoSMIC CLP CL data

Please add the CCLE and CoSMIC CLP Cancer Cell Line (CL) data manually due to legal regulations! Else only the vanilla 62 CellMiner CLs will be available for identification. You can however manually add custom CLs.

'CosmicCLP_MutantExport.tsv.gz', unpack with e.g. gunzip on linux or 7zip on windows from http://cancer.sanger.ac.uk/cell_lines/download

'CCLE_hybrid_capture1650_hg19_NoCommonSNPs_CDS_2012.05.07.tsv' from http://www.broadinstitute.org/ccle/data/browseData

Registration for both websites is without charge and not complicated.

`initiate_canonical_databases( ccle_file = 'path_to_ccle/CCLE_hybrid_capture1650_hg19_NoCommonSNPs_CDS_2012.05.07.tsv', cosmic_file = 'path_to_cosmic/CosmicCLP_MutantExport.tsv.gz' )`

One the initialization succeeds, about 2000 cancer cell line training sample for about 1200 different cancer cell lines are available in the Uniquorn's database. 



Contact: raik.otto@hu-berlin.de

