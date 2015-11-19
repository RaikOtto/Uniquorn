
#' Parse CoSMIC CLP data
parse_cosmic_clp_data = function( parser_path, raw_data  ){
  
  library( "stringr" )
  
  clp_data_path = paste(
    parser_path,
    #'CosmicCLP_CompleteExport.tsv',
    'CosmicCLP_MutantExport.tsv.gz',
    sep = "/"
  )
  
  if ( file.exists( clp_data_path )  ){
    
    message( paste0( "Found CoSMIC file, start parsing :", clp_data_path))
    
    cosmic_coll_class = c( NA, rep("NULL",3), NA, rep("NULL",13), NA, rep("NULL",13) )
    clp_data = read.table( gzfile( clp_data_path ), header = T, fill = T, sep = "\t", colClasses = cosmic_coll_class)
    
    coords = clp_data$Mutation.genome.position
    parse_first_entry = function( cl_string ){ return( as.character( unlist( str_split( cl_string, ":" )  ) )[1] ) }
    chromosomes = as.character( lapply( coords ,FUN = parse_first_entry ) )
    
    parse_second_entry = function( cl_string ){ return( as.character( unlist( str_split( cl_string, "-" )  ) )[2] ) }
    stop_coords = as.character( lapply( coords ,FUN = parse_second_entry ) )
    
    parse_first_entry = function( cl_string ){ return( as.character( unlist( str_split( cl_string, "-" )  ) )[1] ) }
    raw_start_coords  = as.character( lapply( coords ,FUN = parse_first_entry ) )
    parse_second_entry = function( cl_string ){ return( as.character( unlist( str_split( cl_string, ":" )  ) )[2] ) }
    start_coords  = as.character( lapply( raw_start_coords ,FUN = parse_second_entry ) )
    
    new_cosmic_data = data.frame(
      
      "CL_ident" =  as.character( clp_data$Sample.name ),
      "HGNC_symbol" = as.character( clp_data$Gene.name ),
      "Chr" = chromosomes,
      "start" = start_coords,
      "stop" = stop_coords
    )
    
    # filter
    
    #new_cosmic_data = new_cosmic_data[ ! is.na( new_cosmic_data$CL_ident ),]
    new_cosmic_data = new_cosmic_data[ new_cosmic_data$CL_ident != "",]
    new_cosmic_data = new_cosmic_data[ new_cosmic_data$HGNC_symbol != "",]
    new_cosmic_data$CL_ident = paste( new_cosmic_data$CL_ident, "CoSMIC", sep = "_" )
    
    raw_data = rbind(
      raw_data,
      new_cosmic_data
    )
    
    message( paste0( "Parsed CoSMIC file ", clp_data_path ) )
    
  } else {
    
    message( paste0( "Did not find CoSMIC file:", clp_data_path ))
  }
  
  return(raw_data)

}