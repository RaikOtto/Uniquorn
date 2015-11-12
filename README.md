# Younikorn R package

Package to identify cancer cell lines based on their unique somatic mutational fingerprint.

# How to make it work:

1 install R devtools
2 install_github("RaikOtto/Younikorn")
3 install dependencies, in particular "chapmandu2/CancerCellLines"
4 library("Younikorn")
5 Download all required raw data and place them in a folder of you choice, e.g. "Download/raw_data" 
	CCLE www.broadinstitute.org/ccle/data/browseData?conversationPropagation=begin
	Cosmic CLP http://cancer.sanger.ac.uk/cell_lines/download
6 call parse_data_into_younikorn_db( parser_path = "/Download/raw_data" )

Active devlopment

Contact: raik.otto@hu-berlin.de
