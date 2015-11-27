#' Parses the vcf file and predicts the identity of the sample
#' @export
identify_vcf_file = function( vcf_file_path, output_path = "" ){
  
  print( paste0( "Creating fingerprint from VCF file ", vcf_file_path  ) )
  vcf_fingerprint = parse_vcf_file( vcf_file_path )
  
  sim_list_file = paste( system.file("", package = "Uniquorn"), "simlist.RData", sep = "/")
  print( paste0( "Loading similarity data from file ",  sim_list_file )  )
  
  if (! exists("sim_list") )
    attach( sim_list_file  )
  
  print( "Finished loading similarity data. Mapping vcf's fingerprint to all contained fingerprints"  )
  
  path_to_output_db = paste( system.file("", package = "Uniquorn"), "parsed_DB.tab", sep ="/") 
  path_to_output_dict = paste( system.file("", package = "Uniquorn"), "parsed_dict.tab", sep ="/")
  cl_data          = read.table( path_to_output_db,   sep ="\t", header = T )
  all_fingerprints = read.table( path_to_output_dict, sep ="\t", colClasses = c( "character", 'NULL'), header =T )$Fingerprint
  
  mapping = match( all_fingerprints, vcf_fingerprint, nomatch = 0 )
  
  adv = 0
  nr_cls = length( sim_list )
  
  match_fp = function( sim_list_entry ){
    
    stat = round( (adv / as.double(nr_cls)) * 100, 0 )
    adv <<- adv + 1
    
    if ( stat != round( (adv / as.double(nr_cls)) * 100, 0 ) )
      print( paste0( round( (adv / as.double(nr_cls)) * 100, 0 ), " %" ) )
    
    return( sum( as.integer(mapping) & as.integer( sim_list_entry ) ) )
  }
  
  hits = unlist( lapply( sim_list, FUN = match_fp ) )
  candidates = cl_data$CL[ order(hits, decreasing = T)  ]

  
  
  res_tab = data.frame(
    
    "Amount_hits" = hits[order(hits, decreasing = T)],
    "CL_identifier" = candidates
  )
  
  res_tab = res_tab[ res_tab$Amount_hits >= 2  ,]
  
  if ( dim(res_tab)[1] >= 1 ){
    
    print( paste0("Best candidate: ", candidates[1] )  )  
    
    library("stringr")
    
    res_tab[1:5,]
    
    if ( output_path == "" )
      output_path = paste0( vcf_file_path, ".identification.tab" )
    
    print( paste0("Storing information in table: ",output_path ) )
    
    write.table( file = output_path, res_tab, sep ="\t", row.names = F, quote = F  )
    
  } else{
    print( paste0("No CL mutational fingerprint with sufficient similarity found." ) )
  }
  
}