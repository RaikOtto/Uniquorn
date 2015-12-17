#' Loads VCF-based data into the db
#' @export
parse_vcf_file = function( vcf_file_path  ){
  
  suppressPackageStartupMessages(library("stringr"))
  
  if ( file.exists( vcf_file_path ) ){
  
    print( paste0("Reading VCF file: ", vcf_file_path ) )
    
    vcf_handle = read.table( vcf_file_path, sep ="\t", header = F, comment.char = "#" )
    
    vcf_matrix = cbind(
      
      str_replace( str_to_upper( vcf_handle[ , 1 ] ), "CHR", "" ),
      as.integer( vcf_handle[ , 2 ] ),
      str_to_upper( vcf_handle[ , 4 ] ) 
    )
    
    split_add = function( vcf_matrix_row ){
      
      variations = as.character( unlist( str_split( vcf_matrix_row[3], "/" ) ) )
      
      chromosome = rep( vcf_matrix_row[1], length(variations)  )
      start      = as.integer( rep( vcf_matrix_row[2], length(variations)  ) )
      length_vec = unlist( lapply( str_split( variations, "" ), length  ) ) - 1
      end        = start + length_vec
      
      fingerprint = as.character()
      for ( i in 1:length( variations  ) ){ 
        fingerprint = c( fingerprint, paste0( c(chromosome[i], start[i],end[i]), collapse = "_" ) )
      }
      
      return( fingerprint )
    }
    
    fingerprint  = apply( vcf_matrix, FUN = split_add, MARGIN = 1  )
    
    return( fingerprint )
    
  } else {
    
    stop( paste0( "Did not find VCF file: ", vcf_file_path  ) )
  }
}