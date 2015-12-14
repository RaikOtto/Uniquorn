#' Loads VCF-based data into the db
#' @export
parse_vcf_file = function( vcf_file_path  ){
  
  suppressPackageStartupMessages(library("stringr"))
  
  if ( file.exists( vcf_file_path ) ){
  
    message( paste0("Found VCF file: ", vcf_file_path)  )
    
    suppressPackageStartupMessages( require("VariantAnnotation") )
    vcf_handle = readVcf( vcf_file_path, genome = "hg19" )
    ranges     = rowRanges(vcf_handle)@ranges
    start      = as.integer( unclass(ranges)@start )
    len_mut    = as.integer( unclass(ranges)@width )
    end        = start + len_mut - 1
    
    raw_chromosomes = unclass(ranges)@NAMES
    parse_first_entry = function( chromosome_string ){ return( as.character( unlist( str_split( chromosome_string, ":" )  ) )[1] ) }
    chromosomes = as.character( lapply( raw_chromosomes, FUN = parse_first_entry ) )
    chromosomes = str_replace( chromosomes, "chr", "")
    
    vcf_fingerprint_mat = data.frame(
      
      "chromosomes" = chromosomes,
      "start" = start,
      "end" = end
    )
    
    vcf_fingerprint = apply( 
      
      vcf_fingerprint_mat,
      MARGIN = 1,
      FUN = paste,
      collapse = "_",
      sep = ""
    )
    
    vcf_fingerprint = toupper( vcf_fingerprint )
    vcf_fingerprint = as.character( unlist( lapply( vcf_fingerprint, FUN = str_replace_all, " ", "") ) )
    vcf_fingerprint = as.character( unlist( lapply( vcf_fingerprint, FUN = str_replace_all, "CHR", "") ) )
    
    return( vcf_fingerprint )
    
  } else {
    
    stop( paste0( "Did not find VCF file: ", vcf_file_path  ) )
  }
}