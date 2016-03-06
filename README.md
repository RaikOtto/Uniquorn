# Uniquorn R package

Package to identify cancer cell lines (CL)s based on their weighted mutational fingerprint.

# 1 How to make it work:  Quickstart

## Installing the Uniquorn

Start an R session e.g. using RStudio

`source("https://bioconductor.org/biocLite.R")`

`biocLite("Uniquorn")`

## Test run

Here the NCI-60 exome sequenced HT29 Cancer Cell line, reference genome GRCh37/ HG19

`library("Uniquorn")`

`HT29_vcf_file = system.file("extdata/HT29.vcf.gz", package="Uniquorn")`

`ident_result = identify_vcf_file( HT29_vcf_file, ref_gen = "GRCH37"  )`

`head( ident_result )` will show a table with potential identification candidate, how many mutations overall and weighted of the training set have been found and if any training samples have surpased the identification threshold.

### Explanation test data

The HT29 cancer cell line vcf file which contains the somatic mutations of the HT-29 cancer cell line is taken from http://watson.nci.nih.gov/projects/nci60/wes/VCF
The Watson repository contains the same cancer cell line samples as Uniquorn's default training dataset, the original CellMiner panel, available @ http://discover.nci.nih.gov/cellminer/. 
However, the cancer cell line vcf sample was differentely filtered and different algorithms have been used to predict mutations and variations. Therefore, it was used as example set to show the difficulties associated with an identification of cancer cell line sample: Some mutations are found only in the query version (here HT-29 from Watson) and some only in the training dataset (here HT-29 from CellMiner original panel). Moreover, some mutations of the query cancer cell line map to different cancer cell lines in the CellMiner training database what can lead to an incorrect identification of a cancer cell line. 

Therefore, a robust yet sensitive cancer cell line identification algorithm is required.

You will find a file with the ending '_uniquorn_identification.tab' next to the input VCF file if you did not specify the output file path.

# 2 Add CCLE and CoSMIC CLP CL data

Please add the CCLE and CoSMIC CLP Cancer Cell Line (CL) data manually due to legal regulations! Else only the vanilla 62 CellMiner CLs will be available for identification. You can however manually add custom CLs.

'CosmicCLP_MutantExport.tsv.gz', unpack with e.g. gunzip on linux or 7zip on windows from http://cancer.sanger.ac.uk/cell_lines/download

'CCLE_hybrid_capture1650_hg19_NoCommonSNPs_CDS_2012.05.07.tsv' from http://www.broadinstitute.org/ccle/data/browseData

Registration for both websites is without charge and not complicated.

`initiate_canonical_databases( ccle_file = 'path_to_ccle/CCLE_hybrid_capture1650_hg19_NoCommonSNPs_CDS_2012.05.07.tsv', cosmic_file = 'path_to_cosmic/CosmicCLP_MutantExport.tsv.gz', ref_gen = "GRCH37")`

One the initialization succeeds, about 2000 cancer cell line training sample for about 1200 different cancer cell lines are available in the Uniquorn's database.

Note: Currently (January 2016), only the CoSMIC CLP data is available for the reference Genome version GRCh38. It is neccesary, that the reference genome for the training samples is specified if the version is not GRCh37

`initiate_canonical_databases( ccle_file = 'path_to_ccle/CCLE_hybrid_capture1650_hg19_NoCommonSNPs_CDS_2012.05.07.tsv', cosmic_file = 'path_to_cosmic/CosmicCLP_MutantExport.tsv.gz', ref_gen = "GRCH38" )`

# 3 Add training CL samples & utility functions

If you want to identify CL samples not contained in the 'canonical' CL set, you can add your own custom CL samples. These samples will be treated just as the 'canonical' training-datasets from e.g. CCLE. Note however, that it is strongly recommended to add at least 10 sample because overfitting might occur if too little custom training-samples are available. 

`add_custom_vcf_to_database( "path_to_file/my_own_CL_samples.vcf"  )`

Likewisely, if you want to remove the sample:

`remove_custom_vcf_from_database( "Name_of_my_CL_custom_sample"  )`

If you want to see which CLs are contained:

`show_contained_cls( ref_gen = "GRCH37" )`

If you want to know which mutations are overall contained in the training set for a particular database:

`show_contained_mutations( ref_gen = "GRCH37" )`

Same if you want to know which genomic loci are associated with a particular CL:

`show_contained_mutations_for_cl("SF_268_CELLMINER")`

# BED files and Broad Institute IGV visualization

Note as well, that there are BED files for the IGV Browser added as well, so that one can see the 
training, query and missed mutations in the genome. This feature can be switched of by setting the
option `output_bed_file` in the `identify_vcf_file` function `FALSE`.

Contact: raik.otto@hu-berlin.de
