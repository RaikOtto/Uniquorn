
#' Parse CoSMIC CLP data
#' @export
parse_cosmic_clp_data = function( path_to_raw_data  ){
  
  library( "stringr" )
  
  clp_data_path = paste(
    path_to_raw_data,
    'CosmicCLP_CompleteExport.tsv',
    sep = "/"
  )
  
  if ( file.exists( clp_data_path )  ){
    
    clp_data = read.table( clp_data_path, header = T, nrows = 100, fill = T, sep = "\t")
    
    coords = clp_data$Mutation.genome.position
    parse_first_entry = function( cl_string ){ return( as.character( unlist( str_split( cl_string, ":" )  ) )[1] ) }
    chromosomes = as.character( lapply( coords ,FUN = parse_first_entry ) )
    
    parse_second_entry = function( cl_string ){ return( as.character( unlist( str_split( cl_string, "-" )  ) )[2] ) }
    stop_coords = as.character( lapply( coords ,FUN = parse_second_entry ) )
    
    parse_first_entry = function( cl_string ){ return( as.character( unlist( str_split( cl_string, "-" )  ) )[1] ) }
    raw_start_coords  = as.character( lapply( coords ,FUN = parse_first_entry ) )
    parse_second_entry = function( cl_string ){ return( as.character( unlist( str_split( cl_string, ":" )  ) )[2] ) }
    start_coords  = as.character( lapply( raw_start_coords ,FUN = parse_second_entry ) )
    
    raw_data = as.matrix(
      cbind(
        as.character( clp_data$Sample.name ),
        as.character( clp_data$Gene.name ),
        chromosomes,
        start_coords,
        stop_coords
      )
    )
    raw_data = as.data.frame( raw_data )
    colnames( raw_data ) = c( "CL_ident", "HGNC_symbol", "Chr", "start", "stop" )
    
    # filter
    
    raw_data = raw_data[ ! is.na( raw_data$CL_ident ),]
    raw_data = raw_data[ raw_data$CL_ident != "",]
    raw_data = raw_data[ raw_data$HGNC_symbol != "",]
    raw_data$CL_ident = paste( raw_data$CL_ident, "CoSMIC", sep = "_" )
    
    message( paste0( "Parsed file ", hybcappath ) )
    
  } else {
    
    raw_data = matrix( character(), ncol = 5 )
    colnames( raw_data ) = c( "CL_ident", "HGNC_symbol", "Chr", "start", "stop" )
  }
  
  return( raw_data )
}