import argparse, os, gzip

def main_function( parser ):

	cl_db = {}
	cl_dict = {}

	if os.path.exists(parser.cosmic_file):
		
		print 'Parsing cosmic'

		with gzip.open( parser.cosmic_file ) as i_h:

			for i, line in enumerate( i_h ):

				if i == 0: continue

				line = line.strip().split("\t")

				chromosome	= line[ 18 ].split(":")[0].upper().replace("CHR","")
				start_pos	= line[ 18 ].split(":")[1].split("-")[0]
				end_pos		= line[ 18 ].split(":")[1].split("-")[1]
				ident		= line[4].upper()
				ident 		= ident.split("_")[0] + "_COSMIC"
				fingerprint	= "_".join([chromosome,start_pos,end_pos])

				if cl_db.has_key(fingerprint):

					cl_db[fingerprint][ ident ] = True

				else:

					cl_db[ fingerprint ] = {ident:True}
					
				if not cl_dict.has_key(ident):
					
					cl_dict[ident] = {fingerprint:True}
					
				else: cl_dict[ident][fingerprint] = True

	if os.path.exists(parser.ccle_file):
		
		print 'Parsing CCLE'

		with open( parser.ccle_file ) as i_h:

			for i, line in enumerate( i_h ):

				if i == 0: continue

				line = line.strip().split("\t")

				chromosome	= line[4].upper()
				chromosome	= chromosome.replace("CHR","")
				start_pos	= line[5]
				end_pos		= line[6]
				ident		= line[15].upper()
				ident 		= ident.split("_")[0] + "_CCLE"
				fingerprint	= "_".join([chromosome,start_pos,end_pos])

				if cl_db.has_key(fingerprint):

					cl_db[fingerprint][ ident ] = True

				else:
					
					cl_db[ fingerprint ] = {ident:True}
					
				if not cl_dict.has_key(ident):
					
					cl_dict[ident] = {fingerprint:True}
					
				else: cl_dict[ident][fingerprint] = True

	if os.path.exists(parser.cellminer_file):

		print 'Parsing CellMiner'

		cl_names = ["MCF7","MDA_MB_231","HS578T","BT_549","T47D","SF_268","SF_295","SF_539","SNB_19","SNB_75","U251","COLO205","HCC_2998","HCT_116","HCT_15","HT29","KM12","SW_620","CCRF_CEM","HL_60","K_562","MOLT_4","RPMI_8226","SR","LOXIMVI","MALME_3M","M14","SK_MEL_2","SK_MEL_28","SK_MEL_5","UACC_257","UACC_62","MDA_MB_435","MDA_N","A549","EKVX","HOP_62","HOP_92","NCI_H226","NCI_H23","NCI_H322M","NCI_H460","NCI_H522","IGROV1","OVCAR_3","OVCAR_4","OVCAR_5","OVCAR_8","SK_OV_3","NCI_ADR_RES","PC_3","DU_145","786_0","A498","ACHN","CAKI_1","RXF_393","SN12C","TK_10","UO_31"]
		indices = range( len(cl_names ) )

		with open( parser.cellminer_file ) as i_h:

			for i, line in enumerate( i_h ):

				line = line.strip().split("\t")

				if (i <= 10) or (len( line ) < 3 ) : continue

				chromosome	= line[3].upper()
				start_pos	= line[4]
				end_pos		= line[5]
				fingerprint	= "_".join([chromosome,start_pos,end_pos])

				genotype	= line[18:]
				idents = [ ( cl_names[i].upper() + '_CELLMINER' ) for i in indices if genotype[i] != '-' ]

				for ident in idents:

					if cl_db.has_key(fingerprint):

						cl_db[fingerprint][ ident ] = True

					else:

						cl_db[ fingerprint ] = {ident:True}
						
					if not cl_dict.has_key(ident):

						cl_dict[ident] = {fingerprint:True}

					else: cl_dict[ident][fingerprint] = True

	print 'Writing db output'

	with open( parser.output_db, "w" ) as o_h:

		o_h.write( "\t".join( [ "Fingerprint", "CLs" ] ) + "\r\n" )

		for fingerprint in sorted( cl_db.keys() ):

			o_h.write( "\t".join( [ fingerprint, ",".join(cl_db[fingerprint]) ] ) + "\r\n" )

	print 'Writing cl output'

	with open( parser.output_dict, "w" ) as o_h:

		o_h.write( "\t".join( [ "CL", "Fingerprints" ] )  + "\r\n")

		for CL in sorted( cl_dict.keys() ):

			o_h.write( "\t".join( [ CL, ",".join(cl_dict[CL]) ] ) +"\r\n" )

if __name__ == "__main__":

	parser = argparse.ArgumentParser()
	parser.add_argument('-ccle',	'--ccle_file',	type = str, help = 'Input file for ccle',	required = False)
	parser.add_argument('-cosmic',	'--cosmic_file',	type = str, help = 'Input file for cosmic',	required = False)
	parser.add_argument('-cellminer','--cellminer_file',	type = str, help = 'Input  file for cellminer',	required = False)
	parser.add_argument('-o_db',	'--output_db',	type = str, help = 'Path to output_db',	required = True)
	parser.add_argument('-o_dict',	'--output_dict',	type = str, help = 'Path to output dictionary for cl informaiton',	required = True)

	parser = parser.parse_args()
	main_function( parser )
