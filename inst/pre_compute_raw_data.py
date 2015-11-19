def pre_process_matrices( parser ):

    
    
if __name__ == "__main__":

	parser = argparse.ArgumentParser()
	parser.add_argument('-i',   '--input',			type = str,help='Input file', required=False, default="None")
	parser.add_argument('-o',   '--method',			type = str,help='Output file', required=False, default="find")
	parser.add_argument('-type','--program_folder',	type = str, help="folder that the program is located in",required=False,default="./")

	parser = parser.parse_args()
	main_function( parser )
