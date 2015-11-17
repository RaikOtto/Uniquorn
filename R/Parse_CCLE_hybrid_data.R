
#' Parse the CCLE hybrid capture data
parse_ccle_hybrid_data = function( path_to_raw_data, raw_data ){
  
  hybcappath = paste(
    path_to_raw_data,
    'CCLE_hybrid_capture1650_hg19_allVariants_2012.05.07.maf.gz',
    sep = "/"
  )
  
  if ( file.exists( hybcappath )  ){
    
    message( paste0( "Found CCLE file, start parsing :", hybcappath))
    
    data = read.table( gzfile( hybcappath ), header = T, nrows = 1000, fill = T)
    
    parse_first_entry = function( cl_string ){ return( as.character( unlist( str_split( cl_string, "_" )  ) )[1] ) }
    
    new_ccle_data = data.frame(
      
      "CL_ident"    = as.character( lapply( as.character( data$Tumor_Sample_Barcode ), FUN = parse_first_entry ) ),
      "HGNC_symbol" = as.character( data$Hugo_Symbol ),
      "Chr"         = data$Chromosome,
      "start"       = data$Start_position,
      "stop"        = data$End_position
    )
    
    # filter
    
    new_ccle_data = new_ccle_data[ ! is.na( new_ccle_data$CL_ident ),]
    new_ccle_data = new_ccle_data[ new_ccle_data$CL_ident != "",]
    new_ccle_data = new_ccle_data[ new_ccle_data$HGNC_symbol != "",]
    
    new_ccle_data$CL_ident = paste( new_ccle_data$CL_ident, "CCLE", sep = "_" )
    
    raw_data = rbind(
      raw_data,
      new_ccle_data
    )
    
    message( paste0( "Parsed CCLE file ", hybcappath ) )
    
  } else {
    
    message( paste0( "Did not find CCLE file:", hybcappath ))
  }
  
  return(raw_data)

}