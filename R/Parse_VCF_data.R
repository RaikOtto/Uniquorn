#' split_add
#' 
#' @param vcf_matrix_row row of the vcf file
#' @return Transformed entry of vcf file, reduced to start and length
split_add = function( vcf_matrix_row ){
    
    reference  = as.character( unlist( vcf_matrix_row[4] ) )
    length_ref = length( unlist( stringr::str_split( reference  ,"") ))
    variations = as.character( unlist( stringr::str_split( unlist( vcf_matrix_row[5] ), "," ) ) )
    length_variations = length(unlist( stringr::str_split( unlist( vcf_matrix_row[5] ), "," ) ) )
    
    chromosome = rep( as.character( vcf_matrix_row[1] ), length(variations)  )
    start      = as.integer( rep( vcf_matrix_row[2], length(variations)  ) )
    
    fingerprint = as.character()
    for ( i in seq( length_variations ) ){ 
        
        start_var  = as.character( start[i] )
        variation  = variations[i]
        length_var = length( unlist( stringr::str_split( variation,"") ))
        
        length_alt = max( length_var, length_ref )
        end_var    = as.character( as.integer( as.integer(start_var) + length_alt - 1 ) )
        
        chrom      = stringr::str_replace( stringr::str_to_upper( stringr::str_trim( as.character( unlist( chromosome[i] ) ) ) ), "CHR", "" )
        
        fingerprint = c( 
            fingerprint, 
            paste0( c( 
                as.character(chrom),
                as.character(start_var),
                as.character(end_var)
            ),
            collapse = "_" )
        )
    }
    
    return( fingerprint )
}

#' parse_vcf_file
#' 
#' Parses the vcf file and filters all information except for the start and length of variations/ mutations.
#' 
#' @param vcf_file_path Path to the vcf file on the operating system
#' @usage 
#' parse_vcf_file( vcf_file_path )
#' @return Loci-based DNA-mutational fingerprint of the cancer cell line as found in the input VCF file 
parse_vcf_file = function( vcf_file_path  ){
  
  if ( base::file.exists( vcf_file_path ) ){
  
    print( paste0("Reading VCF file: ", vcf_file_path ) )
    
    vcf_handle = utils::read.table( vcf_file_path, sep ="\t", header = FALSE, comment.char = "#", fill = TRUE )

    fingerprint = apply( vcf_handle, FUN = split_add, MARGIN = 1  )
    fingerprint = unique( as.character( unlist(fingerprint ) ) )
    
    return( fingerprint )
    
  } else {
    
    stop( paste0( "Did not find VCF file: ", vcf_file_path  ) )
  }
}