### parse files
#' Parses data into r list variable
#' @export
initiate_uniquorn_database = function( 
    #cosmic_genotype_file = "CosmicCLP_CompleteExport.tsv",
    cosmic_genotype_file = "CosmicCLP_MutantExport.tsv",
    cellminer_genotype_file = 'DNA__Exome_Seq_none.txt',
    ccle_genotype_file = "CCLE_hybrid_capture1650_hg19_allVariants_2012.05.07.maf",
    ucsc_db_snp_file = 'snp142Common.txt',
    ref_gen = "hg19"
  ){
  
  suppressPackageStartupMessages(library("stringr"))
  print( c( "Found CoSMIC: ", file.exists(cosmic_genotype_file) )  )
  print( c( "Found CCLE: ", file.exists(ccle_genotype_file) )  )
  print( c( "Found CellMiner: ", file.exists(cellminer_genotype_file) )  )
  print( c( "Found DbSNP: ", file.exists(ucsc_db_snp_file) )  )
  
  ### pre processing

  #db snp integration
  
  path_to_output_db_unique       = paste( system.file("", package="Uniquorn"), "unique_parsed_DB", sep ="/")
  path_to_output_db_non_unique   = paste( system.file("", package="Uniquorn"), "non_unique_parsed_DB", sep ="/")
  path_to_output_dict_unique     = paste( system.file("", package="Uniquorn"), "unique_parsed_dict", sep ="/")
  path_to_output_dict_non_unique = paste( system.file("", package="Uniquorn"), "non_unique_parsed_dict", sep ="/")

  fingerprint_names_file         = paste0( c( path_to_output_db_unique, "_mut_labels" ), collapse = "")
  fingerprint_names_file_weighted= paste0( c( path_to_output_db_non_unique, "_weighted_mut_labels" ), collapse = "")
  
  stats_file_path          = paste0( c( path_to_output_db_unique,     "_mut_labels_stats" ),    collapse = "" )
  stats_file_path_weighted = paste0( c( path_to_output_db_non_unique, "_weighted_mut_labels_stats" ), collapse = "" )

  path_to_python_dbsnp_python_parser = paste( system.file("", package="Uniquorn"), "parse_db_snp.py", sep ="/")
  path_to_python_dbsnp_python_parser_db = paste( system.file("", package="Uniquorn"), "parse_db_snp_python.pickle", sep ="/")
  path_to_python = paste( system.file("", package="Uniquorn"), "pre_compute_raw_data.py", sep ="/")
  
  if (  file.exists( ucsc_db_snp_file  ) ){
    
    print( paste0( c( "Found DbSNP file",  ucsc_db_snp_file, ", preprocessing."), collapse = " " ) )
    command_line = str_c(
      c(  
        'python', path_to_python_dbsnp_python_parser, 
        "-i", ucsc_db_snp_file,
        "-o", path_to_python_dbsnp_python_parser_db
      ), 
      collapse = " "
    )
      
    if ( ! file.exists(path_to_python_dbsnp_python_parser_db) )
      system( command_line, ignore.stdout = F, intern = F )
    print( "Finished DbSNP pre-processing" )
  }
  
  # unique
  
  command_line = str_c( 
    c(  
      'python',     path_to_python,
      "-ccle ",     ccle_genotype_file,
      "-cosmic ",   cosmic_genotype_file,
      "-cellminer", cellminer_genotype_file,
      "-o_db",      path_to_output_db_unique,
      "-o_dict",    path_to_output_dict_unique,
      "-o_mut_dict",fingerprint_names_file,
      "-i_dbsnp",   path_to_python_dbsnp_python_parser_db,
      "-o_stats_file", stats_file_path,
      "-unique_mode"
    ),
    collapse = " "
  )
  
  if ( ! weighted ){
    system( command_line, ignore.stdout = F, intern = F )
  }

  # non-unique
  
  command_line = str_c( 
    c(  
      'python',        path_to_python,
      "-ccle ",        ccle_genotype_file,
      "-cosmic ",      cosmic_genotype_file,
      "-cellminer",    cellminer_genotype_file,
      "-o_db",         path_to_output_db_non_unique,
      "-o_dict",       path_to_output_dict_non_unique,
      "-o_mut_dict",   fingerprint_names_file_weighted,
      "-o_stats_file", stats_file_path_weighted,
      "-i_dbsnp",      path_to_python_dbsnp_python_parser_db,
      "-filter_frequent"
    ),
    collapse = " "
  )
  
  if ( weighted ){
    system( command_line, ignore.stdout = F, intern = F )
  }
  
  message("Parsing data finished")
  '
  # transform data & load into DB
  
  print( "Loading aggregated fingerprint raw data from all sources"  )
  
  panels = c("CELLMINER","CCLE","COSMIC")
  
  for( panel in panels ){
   
    print(paste0( c(panel) ))
    
    path_to_output_db_panel_unique       = paste0( c( paste( path_to_output_db_unique,      panel, sep ="_" ), ".tab" ), collapse = "")
    path_to_output_db_panel_non_unique   = paste0( c( paste( path_to_output_db_non_unique,  panel, sep ="_" ), ".tab" ), collapse = "")
    
    fingerprint_names_unique_file        = paste0( c( paste( path_to_output_db_unique,      panel, sep ="_" ), "_mut_labels.tab" ), collapse = "")
    fingerprint_names_non_unique_file    = paste0( c( paste( path_to_output_db_non_unique,  panel, sep ="_" ), "_mut_labels.tab" ), collapse = "")
    
    path_to_output_dict_panel_unique     = paste0( c( paste( path_to_output_dict_unique,    panel, sep ="_" ), ".tab" ), collapse = "")
    path_to_output_dict_panel_non_unique = paste0( c( paste( path_to_output_dict_non_unique,panel, sep ="_" ), ".tab" ), collapse = "")
    
    cl_data_unique              = read.table( path_to_output_db_panel_unique,       sep ="\t", header = T )
    cl_data_non_unique          = read.table( path_to_output_db_panel_non_unique,   sep ="\t", header = T )
    
    fingerprint_data_unique     = read.table( path_to_output_dict_panel_unique,     sep ="\t", header = T )
    fingerprint_data_non_unique = read.table( path_to_output_dict_panel_non_unique, sep ="\t", header = T )
    
    sim_list_unique     = create_sim_list( fingerprint_data_unique, cl_data_unique, panel = panel, type = " unique " )
    sim_list_non_unique = create_sim_list( fingerprint_data_non_unique, cl_data_non_unique, panel = panel, type = " non-unique "  )
    
    sim_list_file_unique = paste( system.file("", package = "Uniquorn"), 
      paste0( c("simlist_unique_",panel,".RData"), collapse= "" ), 
      sep = "/"
    )
    
    write.table( x = cbind( as.character( fingerprint_data_unique$Fingerprint ), as.character( fingerprint_data_unique$Weight ) ),
                          fingerprint_names_unique_file, 
                          sep = "\t", quote = F, row.names = F, col.names = F  )
    write.table( x = cbind( as.character( fingerprint_data_non_unique$Fingerprint) , as.character( fingerprint_data_non_unique$Weight )  ), 
                 fingerprint_names_non_unique_file, 
                 sep = "\t", quote = F, row.names = F, col.names = F  )
    
    sim_list_file_non_unique = paste( system.file("", package = "Uniquorn"), 
      paste0( c("simlist_non_unique_",panel,".RData"), 
      collapse= "" ), 
      sep = "/"
    )
    
    print( paste0("Storing similarity information ", sim_list_file_unique)  )
    save( sim_list_unique, file = sim_list_file_unique )
    
    print( paste0("Storing similarity information ", sim_list_file_non_unique)  )
    save( sim_list_non_unique, file = sim_list_file_non_unique )
     
  }
  '
  print("Finished")
}