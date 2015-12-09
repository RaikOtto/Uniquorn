### parse files

#' Parses data into r list variable
#' @export
initiate_uniquorn_database = function( parser_path ){

  library("stringr")
  
  ### pre processing

  clp_data_path = paste(
    
    parser_path,
    'CosmicCLP_CompleteExport.tsv',
    sep = "/"
  )
    
  cellminer_path = paste(
    
    parser_path,
    'DNA__Exome_Seq_none.txt',
    sep = "/"
  )
  
  hybcappath = paste(
    
    parser_path,
    'CCLE_hybrid_capture1650_hg19_NoCommonSNPs_NoNeutralVariants_CDS_2012.05.07.maf',
    sep = "/"
  )
  
  #db snp integration
  
  dbsnp_path = paste( 
    
    parser_path,
    'snp142Common.txt',
    sep = "/"
  )
  
  path_to_output_db = paste( system.file("", package="Uniquorn"), "parsed_DB", sep ="/")  
  path_to_output_dict = paste( system.file("", package="Uniquorn"), "parsed_dict", sep ="/")  
  path_to_python_dbsnp_python_parser = paste( system.file("", package="Uniquorn"), "parse_db_snp.py", sep ="/")
  path_to_python_dbsnp_python_parser_db = paste( system.file("", package="Uniquorn"), "parse_db_snp_python.pickle", sep ="/")
  path_to_python = paste( system.file("", package="Uniquorn"), "pre_compute_raw_data.py", sep ="/")
  
  if (  file.exists( dbsnp_path  ) ){
    
    print( paste0( c( "Found DbSNP file",  dbsnp_path, ", preprocessing."), collapse = " " ) )
    command_line = str_c(
      c(  
        'python', path_to_python_dbsnp_python_parser, 
        "-i", dbsnp_path,
        "-o", path_to_python_dbsnp_python_parser_db
      ), 
      collapse = " "
    )
      
    system( command_line, ignore.stdout = F, intern = F )
    print( "Finished DbSNP pre-processing" )
  }
  
  command_line = str_c( 
    c(  
      'python', path_to_python,
      "-ccle ", hybcappath,
      "-cosmic ", clp_data_path,
      "-cellminer", cellminer_path,
      "-o_db", path_to_output_db,
      "-o_dict", path_to_output_dict,
      "i_dbsnp", path_to_python_dbsnp_python_parser_db
    ),
    collapse = " "
  )
  
  ## aggregate fingerprint with python due to time contrains
  system( command_line, ignore.stdout = F, intern = F )
  
  message("Parsing data finished")
  
  # transform data & load into DB
  
  print( "Loading aggregated fingerprint raw data from all sources"  )
  
  panels = c("CELLMINER","CCLE","COSMIC")
  
  for( panel in panels ){
   
    path_to_output_db_panel   = paste0( c( paste( path_to_output_db,  panel, sep ="_" ), ".tab" ), collapse = "")
    path_to_output_dict_panel = paste0( c( paste( path_to_output_dict,panel, sep ="_" ), ".tab" ), collapse = "")
    
    cl_data          = read.table( path_to_output_db_panel,   sep ="\t", header = T )
    fingerprint_data = read.table( path_to_output_dict_panel, sep ="\t", header = T )
    
    sim_list = create_sim_list( fingerprint_data, cl_data, panel )
    
    sim_list_file = paste( system.file("", package = "Uniquorn"), 
      paste0( c("simlist_",panel,".RData"), collapse= "" ), 
      sep = "/"
    )
    print( paste0("Storing similarity information ", sim_list_file)  )
    
    save( sim_list, file = sim_list_file )
     
  }
  
  print("Finished")
}