### parse files

#' Parses data into the sql database
#' Requires the currently (Nov 2015) github located package "chapmandu2/CancerCellLines"
#' @export
initiate_younikorn_database = function( parser_path, db_path = system.file("", package="Younikorn") ){

  # create db folder
  
  db_path = paste(db_path, "inst",sep="/")
  
  if ( ! dir.exits( db_path )  )
    dir.create( db_path )
  
  db_path  = paste( db_path,"Younikorn.db", sep ="/")
  print( paste0( "Parsing data and storing in db: ",db_path) )

  if ( F ){
  
    # ids data
    
    idspath = system.file(
      "extdata",
      "CellLineIDNormalisationOct15.txt",
      package = "CancerCellLines"
    )
    
    if (file.exists( idspath ) ){
      
      data = read.table( idspath, header = T, sep = "\t")
      dbWriteTable( full_con, "cell_line_ids", data, overwrite = TRUE)
      print( paste0( "Parsed file ", idspath ) )
    }
  
    #infopath
    
    infopath = paste( parser_path, 'CCLE_sample_info_file_2012-10-18.txt', sep = "/" )
    
    if ( file.exists( infopath ) ){
      
      data = read.table( infopath, header = T, sep = "\t")
      colnames( data ) = c(
        "CCLE_name",
        "Primary_cell_name",
        "Cell_line_aliases",
        "Gender",
        "Site_primary",
        "Histology",
        "Hist_subtype1",
        "Notes",
        "Source",
        "Expression_arrays",
        "SNP_arrays",
        "Oncomap",
        "Hybrid_capture_sequencing"
      )
      
      dbWriteTable( 
        full_con,
        "ccle_sampleinfo",
        data,
        overwrite = T
      )
      
      print( paste0( "Parsed file ", infopath ) )
    }
  }
  
  # ccle genotype data

  ccle_raw_data = parse_ccle_hybrid_data( parser_path )

  # Cosmic CLP parsing

  cosmic_raw_data = parse_cosmic_clp_data( parser_path )
  
  # CellMiner NCI60 data
  
  cellminer_raw_data = parse_cellminer_data( parse_path )
  
  message("Parsing data finished")
  
  # create db
  
  require( RSQLite )
  drv = dbDriver("SQLite")
  full_con = dbConnect( drv, dbname = db_path )

}