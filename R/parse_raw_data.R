### parse files

#' Parses data into the sql database
#' Requires the currently (Nov 2015) github located package "chapmandu2/CancerCellLines"
#' @export
initiate_younikorn_database = function( parser_path, db_path = system.file("", package="Younikorn") ){

  library( "stringr" )
  
  raw_data = data.frame(

    "CL_ident" = character(),
    "HGNC_symbol" = character(),
    "Chr" = character(),
    "start" = character(),
    "stop" = character()
  )
  
  # ccle genotype data

  source("./R/Parse_CCLE_hybrid_data.R")
  raw_data = parse_ccle_hybrid_data( parser_path, raw_data )

  # Cosmic CLP parsing

  source("./R/Parse_CoSMIC_data.R")
  raw_data = parse_cosmic_clp_data( parser_path, raw_data )
  
  # CellMiner NCI60 data
  
  source("./R/Parse_nci_60_data.R")
  raw_data = parse_cellminer_data( parser_path, raw_data )
  
  message("Parsing data finished")
  
  # transform data & load into DB
  
  source("./R/Create_similarity_matrix.R")
  similarity_matrix_data = create_similarity_matrix(  )
  
  ### visualization
  
  if (F){
    
    source( "./R/Visualize_matrix.R" )
    visualize_matrix(  )
  }
  
  source("./R/Load_similarity_data_into_db.R")
  load_similarity_data_into_db( similarity_matrix_data, db_path )

}