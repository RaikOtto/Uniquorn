# Younikorn R package

Package to identify cancer cell lines based on their unique somatic mutational fingerprint.

# How to make it work:

1 install('devtools')

2 install_github("RaikOtto/Younikorn")

3 library("Younikorn")

4 install_younikorn_dependencies()

5 Download all required raw data and place them in a folder of you choice, e.g. "Download/raw_data" 
	CCLE www.broadinstitute.org/ccle/data/browseData?conversationPropagation=begin
	Cosmic CLP http://cancer.sanger.ac.uk/cell_lines/download
	CellMiner/ NCI-60 Cls http://discover.nci.nih.gov/cellminerdata/normalizedArchives/nci60_DNA__Exome_Seq_none.zip
	
6 initiate_younikorn_database( parser_path = "/Download/raw_data" )

7 Run toy example

Under active development!

Contact: raik.otto@hu-berlin.de
