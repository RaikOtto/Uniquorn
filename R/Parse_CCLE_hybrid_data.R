
#' Parse the CCLE hybrid capture data
#' @export
parse_ccle_hybrid_data = function( path_to_raw_data  ){
  
  library( "stringr" )
  
  hybcappath = paste(
    path_to_raw_data,
    'CCLE_hybrid_capture1650_hg19_allVariants_2012.05.07.maf.gz',
    sep = "/"
  )
  
  if ( file.exists( hybcappath )  ){
    
    data = read.table( gzfile( hybcappath ), header = T, nrows = 100, fill = T)
    
    parse_first_entry = function( cl_string ){ return( as.character( unlist( str_split( cl_string, "_" )  ) )[1] ) }
    
    raw_data = as.matrix(
      cbind(
        as.character( lapply( as.character( data$Tumor_Sample_Barcode ), FUN = parse_first_entry ) ),
        as.character( data$Hugo_Symbol ),
        data$Chromosome,
        data$Start_position,
        data$End_position
      )
    )
    colnames( raw_data ) = c( "CL_ident", "HGNC_symbol", "Chr", "start", "stop" )
    raw_data = as.data.frame( raw_data )
    # filter
    
    raw_data = raw_data[ ! is.na( raw_data$CL_ident ),]
    raw_data = raw_data[ raw_data$CL_ident != "",]
    raw_data = raw_data[ raw_data$HGNC_symbol != "",]
    
    raw_data$CL_ident = paste( raw_data$CL_ident, "CCLE", sep = "_" )
    
    message( paste0( "Parsed file ", hybcappath ) )
    
  } else {
    
    raw_data = matrix( character(), ncol = 5 )
    colnames( raw_data ) = c( "CL_ident", "HGNC_symbol", "Chr", "start", "stop" )
  }
  
  return( raw_data )
}