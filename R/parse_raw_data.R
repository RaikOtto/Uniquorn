### parse files

#' Parses data into the sql database
#' Requires the currently (Nov 2015) github located package "chapmandu2/CancerCellLines"
#' @export
initiate_younikorn_database = function( parser_path, db_path = system.file("", package="Younikorn") ){

  library("stringr")
  
  ### pre processing
  
  if ( grepl( "/inst", c( db_path )) != T )
    
    db_path = paste( db_path, "inst", sep = "/" )
    
    if ( ! dir.exists( db_path )  )
      dir.create( db_path )
  
    if ( grepl( "Younikorn.db", c( db_path )) != T )
      db_path  = paste( db_path,"Younikorn.db", sep ="/")
  
  clp_data_path = paste(
    
    parser_path,
    #'CosmicCLP_CompleteExport.tsv',
    'CosmicCLP_MutantExport.tsv.gz',
    sep = "/"
  )
    
  cellminer_path = paste(
    
    parser_path,
    'DNA__Exome_Seq_none.txt',
    sep = "/"
  )
  
  hybcappath = paste(
    
    parser_path,
    #'CCLE_hybrid_capture1650_hg19_allVariants_2012.05.07.maf.gz',
    'CCLE_hybrid_capture1650_hg19_NoCommonSNPs_NoNeutralVariants_CDS_2012.05.07.maf',
    sep = "/"
  )
  
  path_to_output_db = paste( system.file("", package="Younikorn"), "inst/parsed_DB.tab", sep ="/")  
  path_to_output_dict = paste( system.file("", package="Younikorn"), "inst/parsed_dict.tab", sep ="/")  
  path_to_python = paste( system.file("", package="Younikorn"), "inst/pre_compute_raw_data.py", sep ="/")
  
  command_line = str_c( 
    c(  
      'python', path_to_python, 
      "-ccle ", hybcappath,
      "-cosmic ", clp_data_path,
      "-cellminer", cellminer_path,
      "-o_db", path_to_output_db,
      "-o_dict", path_to_output_dict
    ), 
  collapse = " " )
  
  ## aggregate fingerprint with python due to time contrains
  system( command_line )
  
  message("Parsing data finished")
  
  # transform data & load into DB
  
  print( "Loading aggregated fingerprint raw data from all sources"  )
  fingerprint_data = read.table( path_to_output_db,   sep ="\t", header = T )
  cl_data          = read.table( path_to_output_dict, sep ="\t", header = T )
  
  source("./R/Create_similarity_matrix.R")
  similarity_matrix_data = create_similarity_matrix( fingerprint_data, cl_data )
  
  ### visualization
  
  if (F){
    
    source( "./R/Visualize_matrix.R" )
    visualize_matrix(  )
  }
  
  source("./R/Load_similarity_data_into_db.R")
  load_similarity_data_into_db( similarity_matrix_data, db_path )

}