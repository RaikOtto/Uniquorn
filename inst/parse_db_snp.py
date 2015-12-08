"""

UCSC dbsnp files, such as snp142Common.txt

"""
epilog="""Raik Otto <raik.otto@hu-berlin.de> 20.1.2015"""

import argparse, os, cPickle as pickle

def parse_dbsnp( parser ):
	
	print "Commenced parsing "+parser.db_snp_file
	
	db_snp_d = {}
	verbose_d = {}
	
	with open( parser.db_snp_file ) as i_h:
		
		for line in i_h:
			
			line = line.split("\t")
			
			chrom = line[1].replace("chr","")
			if "_" in chrom: continue
			
			if not verbose_d.has_key(chrom):
				
				verbose_d[chrom] = True
				print chrom
			
			start = line[2]
			end   = line[3]
			signature = "_".join( [ chrom, start, end ] )

			db_snp_d[signature ] = True
			
	print "Dumping dat to " + parser.pickle_output_file
	pickler = pickle.Pickler( file( "db_snp_python_obj.pickle" , "wb"), -1)
	pickler.dump( db_snp_d )
	print "Dumping finished"

if __name__ == "__main__":

	parser = argparse.ArgumentParser()
	parser.add_argument('-i',	'--db_snp_file',		type = str, help = 'Input file for db snp, tab separated UCSC',	required = True)
	parser.add_argument('-o',	'--pickle_output_file',	type = str, help = 'pickle_output_file',	required = False, default = "./db_snp_python_obj.pickle")

	parser = parser.parse_args()
	parse_dbsnp( parser )

