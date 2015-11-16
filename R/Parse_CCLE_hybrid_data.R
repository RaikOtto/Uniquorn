
#' Parse the CCLE hybrid capture data
#' @export
parse_ccle_hybrid_data = function( path_to_raw_data  ){
  
  hybcappath = paste( 
    path_to_raw_data,
    'CCLE_hybrid_capture1650_hg19_NoCommonSNPs_NoNeutralVariants_CDS_2012.05.07.xlsx',
    sep = "/"
  )
  
  if ( file.exists( hybcappath ) ){
    
    require( readxl )
    
    data = read_excel(
      
      hybcappath,
      col_types = c( 
        "text",
        "numeric",
        "text",
        "numeric",
        "text",
        "numeric",
        "numeric",
        rep("text",44)
      )
    )
    
    print(
      paste0(
        "Parsed file ",
        hybcappath
      )
    )
  }
  
  return( "data_matrix" )
}