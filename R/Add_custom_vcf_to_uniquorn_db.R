#' Add a custum vcf file to the three existing cancer cell line panels
#' @export
add_custom_vcf_to_uniquorn_db = function( vcf_file_path, ref_gen = "hg19" ){
    
  library("stringr")
  panels = c("CELLMINER","CCLE","COSMIC"), 
  
  print( paste0( "Creating fingerprint from custom VCF file ", vcf_file_path  ) )
  vcf_fingerprint = parse_vcf_file( vcf_file_path )
    
  sim_list_store_mat = matrix( list(), ncol = length( panels ) )
  colnames(sim_list_store_mat) = panels
  sim_list_stats_store_mat = matrix( list(), ncol = length( panels ) )
  colnames(sim_list_stats_store_mat) = panels
  
}