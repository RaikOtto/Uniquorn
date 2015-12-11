"""

Input: The three mutation/ variantion holding text files from CCLE, CoSMIC CLP and CellMiner NCI-60 
Outout: Dataset holding unique fingerprints. Uniqueness only valid within a particular dataset

"""
epilog="""Raik Otto <raik.otto@hu-berlin.de> 20.1.2015"""

import argparse, os, cPickle as pickle, operator, math, functools

def percentile( N, percent, key = lambda x : x ):
	"""
	## {{{ http://code.activestate.com/recipes/511478/ (r1)
	Find the percentile of a list of values.
	
	@parameter N - is a list of values. Note N MUST BE already sorted.
	@parameter percent - a float value from 0.0 to 1.0.
	@parameter key - optional key function to compute value from each element of N.
	
	@return - the percentile of the values
	
	Authorship: Wai Yip Tung
	"""
	
	if not N: return None
	
	k = (len(N)-1) * percent
	f = math.floor(k)
	c = math.ceil(k)
	
	if f == c: return key(N[int(k)])
	
	d0 = key(N[int(f)]) * (c-k)
	d1 = key(N[int(c)]) * (k-f)
	
	return d0+d1

def load_data( parser ):

	cl_db   = { 'CCLE':{}, 'COSMIC':{}, 'CellMiner':{} } # stores cell lines for every mutation
	cl_dict = { 'CCLE':{}, 'COSMIC':{}, 'CellMiner':{} } # stores mutation for every cell line
	stat_d  = { 'CCLE':{}, 'COSMIC':{}, 'CellMiner':{} } # counts how often you see a mutation to filter for the often occuring ones
	seen_d  = { 'CCLE':{}, 'COSMIC':{}, 'CellMiner':{} } # counts how often you see a mutation 

	cellminer_cl_names = ["MCF7","MDA_MB_231","HS578T","BT_549","T47D","SF_268","SF_295","SF_539","SNB_19","SNB_75","U251","COLO205","HCC_2998","HCT_116","HCT_15","HT29","KM12","SW_620","CCRF_CEM","HL_60","K_562","MOLT_4","RPMI_8226","SR","LOXIMVI","MALME_3M","M14","SK_MEL_2","SK_MEL_28","SK_MEL_5","UACC_257","UACC_62","MDA_MB_435","MDA_N","A549","EKVX","HOP_62","HOP_92","NCI_H226","NCI_H23","NCI_H322M","NCI_H460","NCI_H522","IGROV1","OVCAR_3","OVCAR_4","OVCAR_5","OVCAR_8","SK_OV_3","NCI_ADR_RES","PC_3","DU_145","786_0","A498","ACHN","CAKI_1","RXF_393","SN12C","TK_10","UO_31"]
	cellminer_indices = range( len( cellminer_cl_names ) )

	# pre loading dbsnp

	db_snp_mode = os.path.isfile( parser.pickle_dbsnp_file )

	if db_snp_mode:

		print "Reading DbSNP dump " + parser.pickle_dbsnp_file
		unpickler = pickle.Unpickler( file( parser.pickle_dbsnp_file , "r"))
		db_snp_d = unpickler.load()
		snp_count = 0
		print "Finished unpickling"
		
	else:
		
		print "Did not find UCSC python DbSNP file, e.g. "

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

					# query for membership - non_unique version with db SNP

					if db_snp_mode: # check if db snp dictionary is available

						if db_snp_d.has_key(fingerprint): 

							snp_count = snp_count + 1
							continue

					for cl_ident in ident_list:

						if parser.unique_mode:

							if ( not seen_d[ type_panel ].has_key( fingerprint ) ):

								seen_d[ type_panel ][fingerprint] = True
								if not cl_db[   type_panel ].has_key( cl_ident ):     cl_db[   type_panel ][ cl_ident ] = {}
								if not cl_dict[   type_panel ].has_key( fingerprint ):cl_dict[ type_panel ][ fingerprint ] = {}
								if not stat_d[   type_panel ].has_key( fingerprint ): stat_d[  type_panel ][ fingerprint ] = 0

								cl_dict[ type_panel ][ fingerprint ][ cl_ident ] = True
								cl_db[   type_panel ][ cl_ident    ][ fingerprint ] = True # save fingerprint for cl
								stat_d[ type_panel ] [ fingerprint ] = stat_d[ type_panel ] [ fingerprint ] + 1

							else:


								if cl_db[ type_panel ].has_key(cl_ident):
									if cl_db[   type_panel ][ cl_ident ].has_key(fingerprint):

										del cl_db[ type_panel ][ cl_ident ][ fingerprint ]


								if cl_dict[   type_panel ].has_key(fingerprint):

									del cl_dict[ type_panel ][ fingerprint ]

								if parser.unique_mode and stat_d[ type_panel ].has_key(fingerprint) : del stat_d[ type_panel ] [ fingerprint ]
								
						else: 

							if not cl_db[   type_panel ].has_key( cl_ident ):     cl_db[   type_panel ][ cl_ident ] = {}
							if not cl_dict[   type_panel ].has_key( fingerprint ):cl_dict[ type_panel ][ fingerprint ] = {}
							if not stat_d[   type_panel ].has_key( fingerprint ): stat_d[  type_panel ][ fingerprint ] = 0

							cl_dict[ type_panel ][ fingerprint ][ cl_ident ] = True
							cl_db[   type_panel ][ cl_ident    ][ fingerprint ] = True # save fingerprint for cl
							stat_d[ type_panel ] [ fingerprint ] = stat_d[ type_panel ] [ fingerprint ] + 1

		# filter part
		"""
		print "Before filtering ", type_panel, ": " , len(stat_d[type_panel].keys())
		
		if ( not parser.retain_frequent_mutations ):

			print( 'Filtering the most frequent mutations' )
			
			perc = .1
			cutoff = percentile( N = stat_d[ type_panel ].values()  , percent = perc )

			print round(cutoff,2), " percentile filter: " , perc, type_panel

			for fingerprint in cl_dict[ type_panel ].keys():

				if float(stat_d[ type_panel ][fingerprint]) > float(cutoff):

					del cl_dict[type_panel][fingerprint]
					del stat_d[type_panel][fingerprint]

			print "After filtering ", type_panel, ": " , len(stat_d[type_panel].keys())
		"""

	if db_snp_mode: print "Excluded " + str(snp_count) + " many SNPs"

	print( 'Writing output' )

	for type_panel in [ "CCLE", "CellMiner", "COSMIC" ]:

		if type_panel == 'COSMIC':

			in_file = parser.cosmic_file

		elif type_panel == 'CCLE':

			in_file = parser.ccle_file

		elif type_panel == 'CellMiner':

			in_file = parser.cellminer_file

		if os.path.exists( in_file ):


			with open( parser.output_mut_dict + "_" + type_panel + ".tab", "w" ) as o_h:

				o_h.write( "\t".join( [ "Fingerprint", "CL", "Count" ] ) + "\r\n" )

				for fingerprint in sorted( cl_dict[ type_panel ].keys() ):

					member_cl = cl_dict[ type_panel ][fingerprint].keys()[0]
					member_count = str( len( cl_db[   type_panel ][ member_cl ].keys() ) )
					o_h.write( "\t".join( [ fingerprint, member_cl, member_count ] ) + "\r\n" )

			"""
			with open( parser.output_dict + "_" + type_panel + ".tab", "w" ) as o_h:

				o_h.write( "\t".join( [ "Fingerprint", "CLs", "Weight" ] ) + "\r\n" )

				for fingerprint in sorted( cl_dict[ type_panel ].keys() ):

					member_cls = ",".join( cl_dict[ type_panel ][fingerprint].keys() )
					o_h.write( "\t".join( [ fingerprint, member_cls, str( round( 1.0 / len(cl_dict[ type_panel ][fingerprint].keys()), 3 ) ) ] ) + "\r\n" )


			with open( parser.output_db + "_" + type_panel + ".tab", "w" ) as o_h:

				o_h.write( "\t".join( [ "CL", "Fingerprints" ] )  + "\r\n")


				for CL in sorted( cl_db[ type_panel ].keys() ):

					found_mutations = [ mut_ident for mut_ident in cl_db[ type_panel ][ CL ] if cl_dict[type_panel].has_key(mut_ident) ]
					o_h.write( "\t".join( [ CL, ",".join( found_mutations ) ] ) +"\r\n" )
			"""

	print("Finished data parsing")

if __name__ == "__main__":

	parser = argparse.ArgumentParser()
	parser.add_argument('-ccle',	'--ccle_file',		type = str, help = 'Input file for ccle',				required = False)
	parser.add_argument('-cosmic',	'--cosmic_file',	type = str, help = 'Input file for cosmic',		required = False)
	parser.add_argument('-cellminer','--cellminer_file',type = str, help = 'Input  file for cellminer',	required = False)
	#parser.add_argument('-db',	'--db_path',	type = str, help = 'Path to output_db',	required = True)
	parser.add_argument('-i_dbsnp',	'--pickle_dbsnp_file',type = str, help = 'pickle_output_file of db snp',	required = False, default = "")
	parser.add_argument('-o_db',	'--output_db',		type = str, help = 'Path to output_db',	required = True)
	parser.add_argument('-o_dict',	'--output_dict',	type = str, help = 'Path to output dictionary for cl information',	required = True)
	parser.add_argument('-o_mut_dict',	'--output_mut_dict',	type = str, help = 'Path to output dictionary for cl information',	required = True)
	parser.add_argument('-retain_frequent',	'--retain_frequent_mutations', action='store_true',	help = 'Filter 50% most frequent mutations' )
	parser.add_argument('-unique_mode',	'--unique_mode', action='store_true',	help = 'Only use unique mutations', default = False )

	parser = parser.parse_args()
	load_data( parser )
