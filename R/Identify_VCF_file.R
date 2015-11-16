#' Parses the vcf file and predicts the identity of the sample
#' @export
identify_vcf_file = function( vcf_file_path ){
  
  require( "VariantAnnotation" )
  parse_vcf_file
    
}