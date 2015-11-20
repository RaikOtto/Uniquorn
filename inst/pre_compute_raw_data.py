from __future__ import print_function
"""

Input: The three mutation/ variantion holding text files from CCLE, CoSMIC CLP and CellMiner NCI-60 
Outout: Dataset holding unique fingerprints. Uniqueness only valid within a particular dataset

"""
epilog="""Raik Otto <raik.otto@hu-berlin.de> 20.1.2015"""

import argparse, os, sqlite3


def load_data( parser ):

	cl_db   = { 'CCLE':{}, 'COSMIC':{}, 'CELLMINER':{} } # stores mutation for every cell line
	cl_dict = { 'CCLE':{}, 'COSMIC':{}, 'CELLMINER':{} } # stores cell lines for every mutation
	seen_db = { 'CCLE':{}, 'COSMIC':{}, 'CELLMINER':{} } # which mutations have already been found

	cellminer_cl_names = ["MCF7","MDA_MB_231","HS578T","BT_549","T47D","SF_268","SF_295","SF_539","SNB_19","SNB_75","U251","COLO205","HCC_2998","HCT_116","HCT_15","HT29","KM12","SW_620","CCRF_CEM","HL_60","K_562","MOLT_4","RPMI_8226","SR","LOXIMVI","MALME_3M","M14","SK_MEL_2","SK_MEL_28","SK_MEL_5","UACC_257","UACC_62","MDA_MB_435","MDA_N","A549","EKVX","HOP_62","HOP_92","NCI_H226","NCI_H23","NCI_H322M","NCI_H460","NCI_H522","IGROV1","OVCAR_3","OVCAR_4","OVCAR_5","OVCAR_8","SK_OV_3","NCI_ADR_RES","PC_3","DU_145","786_0","A498","ACHN","CAKI_1","RXF_393","SN12C","TK_10","UO_31"]
	cellminer_indices = range( len( cellminer_cl_names ) )

	for type_panel in [ "CCLE","COSMIC", "CELLMINER" ]:

		if type_panel == 'COSMIC':

			in_file = parser.cosmic_file

		elif type_panel == 'CCLE':

			in_file = parser.ccle_file

		elif type_panel == 'CELLMINER':

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
						ident_list 	= [ ( cl_names[i].upper() + '_CELLMINER' ) for i in indices if genotype[i] != '-' ]

					fingerprint	= "_".join([chromosome,start_pos,end_pos]) # definition fp

					# query for membership - unique version

					for cl_ident in ident_list:

						if not seen_db[ type_panel ].has_key( fingerprint ):

							if not cl_db[   type_panel ].has_key( cl_ident ): cl_db[   type_panel ][ cl_ident ] = {}

							seen_db[ type_panel ][ fingerprint ] = cl_ident # store as seen
							cl_db[   type_panel ][ cl_ident    ][ fingerprint ] = True # save fingerprint for cl
							cl_dict[ type_panel ][ fingerprint ] = cl_ident

						else:

							seen_cl_ident = seen_db[ type_panel ][ fingerprint ] # find cl for which the fingerprint has been observed

							if cl_db[   type_panel ][ seen_cl_ident ].has_key(fingerprint):

								del cl_db[   type_panel ][ seen_cl_ident ][ fingerprint ]


							if cl_dict[   type_panel ].has_key(fingerprint):

								del cl_dict[ type_panel ][ fingerprint ]

	"""
	# fp information

	conn = sqlite3.connect( parser.db_path)
	c = conn.cursor()
	c.execute('CREATE TABLE IF NOT EXISTS fingerprint_information_table (fingerprint text, CLs text)')

	insertion = []

	for type_panel in [ "CCLE", "CELLMINER", "COSMIC" ]:

		print ("Loading to DB fingerprint information of: "+type_panel)

		for fp in sorted( cl_dict[ type_panel ].keys() ):

			cls = cl_dict[ type_panel ][ fp]
			insertion = insertion + [ (fp, cls) ]

	c.executemany('INSERT INTO fingerprint_information_table VALUES (?,?)', insertion)

	# cl information

	c.execute('CREATE TABLE IF NOT EXISTS cell_line_information_table (CL text, fingerprints text)')

	insertion = []

	for type_panel in [ "CCLE", "CELLMINER", "COSMIC" ]:

		print ("Loading to DB cl information of: "+type_panel)

		for cl in sorted( cl_db[ type_panel ].keys() ):

			fps = ",".join( cl_db[ type_panel ][ cl ] )
			insertion = insertion + [ (str(cl), str(fps) ) ]

	c.executemany('INSERT INTO cell_line_information_table VALUES (?,?)', insertion)
	conn.commit()

	conn.close()


	"""

	print( 'Writing db output' )

	with open( parser.output_db, "w" ) as o_h:

		o_h.write( "\t".join( [ "Fingerprint", "CLs" ] ) + "\r\n" )

		for type_panel in [ "CCLE", "CELLMINER", "COSMIC" ]:

			for fingerprint in sorted( cl_db[ type_panel ].keys() ):

				#o_h.write( "\t".join( [ fingerprint, ",".join( cl_db[ type_panel ][fingerprint] ) ] ) + "\r\n" )
				o_h.write( "\t".join( [ fingerprint, cl_db[ type_panel ][fingerprint] ] ) + "\r\n" )

	print( 'Writing cl output' )

	with open( parser.output_dict, "w" ) as o_h:

		o_h.write( "\t".join( [ "CL", "Fingerprints" ] )  + "\r\n")

		for type_panel in [ "CCLE", "CELLMINER", "COSMIC" ]:

			for CL in sorted( cl_dict[ type_panel ].keys() ):

				#o_h.write( "\t".join( [ CL, ",".join( cl_dict[ type_panel ][ CL ]) ] ) +"\r\n" )
				o_h.write( "\t".join( [ CL, ",".join( cl_dict[ type_panel ][ CL ]) ] ) +"\r\n" )

	print("Finished data parsing")

if __name__ == "__main__":

	parser = argparse.ArgumentParser()
	parser.add_argument('-ccle',	'--ccle_file',	type = str, help = 'Input file for ccle',	required = False)
	parser.add_argument('-cosmic',	'--cosmic_file',	type = str, help = 'Input file for cosmic',	required = False)
	parser.add_argument('-cellminer','--cellminer_file',	type = str, help = 'Input  file for cellminer',	required = False)
	#parser.add_argument('-db',	'--db_path',	type = str, help = 'Path to output_db',	required = True)
	parser.add_argument('-o_db',	'--output_db',	type = str, help = 'Path to output_db',	required = True)
	parser.add_argument('-o_dict',	'--output_dict',	type = str, help = 'Path to output dictionary for cl informaiton',	required = True)

	parser = parser.parse_args()
	load_data( parser )
