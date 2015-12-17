#' Adds a custom vcf file to the three existing cancer cell line panels
#' @export
add_custom_vcf_to_uniquorn_db = function( vcf_file_path, ref_gen = "hg19", name_vcf = "" ){
    
  library("stringr")
  panels = c("CELLMINER","CCLE","COSMIC")
  ref_gen_path = paste( system.file("", package = "Uniquorn"), ref_gen, sep ="/" )
  
  if ( file.exists(vcf_file_path)  ){
  
    print( paste0( "Creating fingerprint from custom VCF file ", vcf_file_path  ) )
    
    if (name_vcf == ""){
      
      vcf_identifier = as.character( tail( unlist( str_split( vcf_file_path, "/" ) ), 1) )
      vcf_identifier = head( unlist( str_split( vcf_identifier, ".vcf|.VCF" ) )  , 1)
      #vcf_identifier = 
      
    } else {
      
      vcf_identifier = name_vcf
    }
    
    vcf_fingerprint = parse_vcf_file( vcf_file_path )
    message( paste0( "Found file ",  vcf_file_path )  )
    
    for( panel in panels ){
      
      message( panel  )
      
      sim_list_file       = paste0( c( ref_gen_path, "/", "Fingerprint_",       panel, ".tab" ), collapse = "" )
      sim_list_stats_file = paste0( c( ref_gen_path, "/", "Fingerprint_stats_", panel, ".tab" ), collapse = "" )
      
      sim_list       = read.table( sim_list_file, sep = "\t", header = T)
      sim_list_stats = read.table( sim_list_stats_file, sep = "\t", header = T)
      
      #if ( !   )
      
      list_of_cls       = sort( unique( sim_list$CL ) )
      nr_cls            = length( list_of_cls  ) # amount cls
      found_mut_mapping = which( sim_list$Fingerprint %in% vcf_fingerprint ) # mapping
      cl_match_stats    = match( list_of_cls, sim_list_stats$CL   ) # mapping
    }
  } else {
    
    message( paste0( "Did not find file: ", vcf_file_path)  )
  } 
}