#' Parses the vcf file and predicts the identity of the sample
#' @export
identify_vcf_file = function( vcf_file_path ){
  
  vcf_fingerprint = parse_vcf_file( vcf_file_path )
  
  identify_vcf_fingerprint( vcf_fingerprint  )
    
}