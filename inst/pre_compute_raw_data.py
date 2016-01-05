"""

Input: The three mutation/ variantion holding text files from CCLE, CoSMIC CLP and CellMiner NCI-60 
Outout: Dataset holding unique fingerprints. Uniqueness only valid within a particular dataset

"""
epilog="""Raik Otto <raik.otto@hu-berlin.de> 20.1.2015"""

import argparse, os, cPickle as pickle, operator, math, functools

def load_data( parser ):

	cl_db   = { 'CCLE':{}, 'COSMIC':{}, 'CellMiner':{} } # stores cell lines for every mutation
	cl_dict = { 'CCLE':{}, 'COSMIC':{}, 'CellMiner':{} } # stores mutation for every cell line

	cellminer_cl_names = ["MCF7","MDA_MB_231","HS578T","BT_549","T47D","SF_268","SF_295","SF_539","SNB_19","SNB_75","U251","COLO205","HCC_2998","HCT_116","HCT_15","HT29","KM12","SW_620","CCRF_CEM","HL_60","K_562","MOLT_4","RPMI_8226","SR","LOXIMVI","MALME_3M","M14","SK_MEL_2","SK_MEL_28","SK_MEL_5","UACC_257","UACC_62","MDA_MB_435","MDA_N","A549","EKVX","HOP_62","HOP_92","NCI_H226","NCI_H23","NCI_H322M","NCI_H460","NCI_H522","IGROV1","OVCAR_3","OVCAR_4","OVCAR_5","OVCAR_8","SK_OV_3","NCI_ADR_RES","PC_3","DU_145","786_0","A498","ACHN","CAKI_1","RXF_393","SN12C","TK_10","UO_31"]
	cellminer_indices = range( len( cellminer_cl_names ) )

	# pre loading dbsnp

	for type_panel in [ "CellMiner", "CCLE", "COSMIC" ]:

		if type_panel == 'COSMIC':

			in_file = parser.cosmic_file

		elif type_panel == 'CCLE':

			in_file = parser.ccle_file

		elif type_panel == 'CellMiner':

			in_file = parser.cellminer_file

		# iterate over mutations

		if os.path.exists( in_file ):

			print( 'Parsing '+ in_file )

			with open( in_file ) as i_h:

				for i, line in enumerate( i_h ):

					if i == 0: continue

					line = line.strip().split("\t")

					if type_panel == 'COSMIC':

						if ( (i == 0) or ( len( line ) < 3 ) or ( len( line[ 23 ].split(":") ) < 2 ) ) : continue

						#print( line)
						chromosome	= line[ 23 ].split(":")[0].upper().replace("CHR","")
						start_pos	= line[ 23 ].split(":")[1].split("-")[0]
						end_pos		= line[ 23 ].split(":")[1].split("-")[1]
						cl_ident	= line[4].upper()
						cl_ident	= cl_ident.split("_")[0] + "_COSMIC"
						ident_list  = [ cl_ident ]

					elif type_panel == 'CCLE':

						if (i == 0) or (len( line ) < 3 ) : continue

						chromosome	= line[4].upper()
						chromosome	= chromosome.replace("CHR","")
						start_pos	= line[5]
						end_pos		= line[6]
						cl_ident	= line[15].upper()
						cl_ident 	= cl_ident.split("_")[0] + "_CCLE"
						ident_list  = [ cl_ident ]

					elif type_panel == 'CellMiner':

						if (i <= 10) or (len( line ) < 3 ) : continue

						chromosome	= line[3].upper()
						start_pos	= line[4]
						end_pos		= line[5]

						genotype	= line[18:]
						ident_list 	= [ ( cellminer_cl_names[i].upper() + '_CELLMINER' ) for i in cellminer_indices if genotype[i] != '-' ]

					fingerprint	= "_".join([chromosome,start_pos,end_pos]) # definition fp

					for cl_ident in ident_list:

						if not cl_db[   type_panel ].has_key( cl_ident ):     cl_db[   type_panel ][ cl_ident ] = {}
						if not cl_dict[   type_panel ].has_key( fingerprint ):cl_dict[ type_panel ][ fingerprint ] = {}

						cl_dict[ type_panel ][ fingerprint ][ cl_ident ] = True
						cl_db[   type_panel ][ cl_ident    ][ fingerprint ] = True # save fingerprint for cl

			# filter part

			print( 'Writing output' )

			with open( parser.o_db_path + "/Fingerprint_" +  type_panel + ".tab", "w" ) as o_h:

					o_h.write( "\t".join( [ "Fingerprint", "CL" ] ) + "\r\n" )

					for fingerprint in sorted( cl_dict[ type_panel ].keys() ):

						member_cls = cl_dict[ type_panel ][ fingerprint ].keys()

						for member_cl in member_cls:

							o_h.write( "\t".join( [ fingerprint, member_cl ] ) + "\r\n" )

		else:

			print "Could not find file: " + in_file

	print("Finished pre-calculation of raw data, parsing pre-calculated data")

if __name__ == "__main__":

	parser = argparse.ArgumentParser()
	parser.add_argument('-ccle',	'--ccle_file',				type = str, help = 'Input file for ccle',			required = False)
	parser.add_argument('-cosmic',	'--cosmic_file',			type = str, help = 'Input file for cosmic',			required = False)
	parser.add_argument('-cellminer','--cellminer_file',		type = str, help = 'Input  file for cellminer',		required = False)
	parser.add_argument('-o_db_path',	'--o_db_path',	type = str, help = 'Path to output_db',				required = True)

	parser = parser.parse_args()
	load_data( parser )
