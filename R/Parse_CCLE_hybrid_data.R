
#' Parse the CCLE hybrid capture data
parse_ccle_hybrid_data = function( parser_path, raw_data ){
  
  hybcappath = paste(
    parser_path,
    #'CCLE_hybrid_capture1650_hg19_allVariants_2012.05.07.maf.gz',
    'CCLE_hybrid_capture1650_hg19_NoCommonSNPs_NoNeutralVariants_CDS_2012.05.07.maf',
    sep = "/"
  )
  
  if ( file.exists( hybcappath )  ){
    
    message( paste0( "Found CCLE file, start parsing :", hybcappath))
    
    ccle_col_clases = c( NA, rep("NULL",3), NA, NA,NA, rep("NULL",8), NA, rep("NULL",35) )

    hybrid_data = read.table(
      check.names = F,
      na.strings=c(""),
      hybcappath,
      header = T,
      fill = T,
      #sep = "\t",
      colClasses = ccle_col_clases,
      nrows = 10**3
    )
    
    parse_first_entry = function( cl_string ){ return( as.character( unlist( str_split( cl_string, "_" )  ) )[1] ) }
    
    new_ccle_data = data.frame(
      
      "CL_ident"    = as.character( lapply( as.character( hybrid_data$Tumor_Sample_Barcode ), FUN = parse_first_entry ) ),
      "HGNC_symbol" = as.character( hybrid_data$Hugo_Symbol ),
      "Chr"         = hybrid_data$Chromosome,
      "start"       = hybrid_data$Start_position,
      "stop"        = hybrid_data$End_position
    )
    
    # filter
    
    new_ccle_data = new_ccle_data[ ! is.na( new_ccle_data$CL_ident ),]
    new_ccle_data = new_ccle_data[ new_ccle_data$CL_ident != "null"  ,]
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