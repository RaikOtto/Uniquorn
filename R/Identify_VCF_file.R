#' Parses the vcf file and predicts the identity of the sample
#' @export
identify_vcf_file = function( vcf_file_path, output_path = "" ){
  
  library("stringr")
  
  print( paste0( "Creating fingerprint from VCF file ", vcf_file_path  ) )
  vcf_fingerprint = parse_vcf_file( vcf_file_path )
  
  panels = c("CELLMINER","CCLE","COSMIC")
  types = c( "unique", "non_unique")
  
  for( type in types  ){
  for( panel in panels ){
  
    if ( output_path == ""  ){
      output_path_panel = paste( vcf_file_path, paste0( c( paste( panel, type, sep ="_" ),".tab"), collapse = ""), sep ="_ident_")
    } else {
      output_path_panel = paste( paste0( output_path, tail( unlist(str_split(vcf_file_path, "/" ) ), 1 )  ), paste0( c( paste( panel, type, sep ="_" ),".tab"), collapse = ""), sep ="_ident_")
    }
    
    
    sim_list_file = paste( system.file("", package = "Uniquorn"), 
       paste0( 
         c(
           paste0( c("simlist", type,panel), collapse = "_" ),
           ".RData"
         ), 
         collapse= ""
       ), 
       sep = "/"
    )

    fingerprint_names_file = paste( system.file("", package = "Uniquorn"), 
        paste0( 
           c(
             type,
             "_parsed_DB_",
             panel,
             "_mut_labels.tab"
            ), collapse= ""
        ), 
        sep = "/"
    )
    
    print( paste0( "Loading similarity data from file ",  sim_list_file )  )
    
    #if (! exists("sim_list") )
    attach( sim_list_file  )
    
    if (type == "unique"){
      sim_list = sim_list_unique
    } else {
      sim_list = sim_list_non_unique
    }
    all_fingerprints_in_cl_set = read.table( fingerprint_names_file, sep ="\t", header= F )[,1]
    
    print( "Finished loading similarity data. Mapping vcf's fingerprint to all contained fingerprints"  )
    
    mapping = match( all_fingerprints_in_cl_set, vcf_fingerprint, nomatch = 0 )
    
    adv = 0
    nr_cls = length( sim_list )
    
    match_fp = function( sim_list_entry ){
      
      stat = round( (adv / as.double(nr_cls)) * 100, 0 )
      adv <<- adv + 1
      
      if ( stat != round( (adv / as.double(nr_cls)) * 100, 0 ) )
        print( paste0( c( panel,  round( (adv / as.double(nr_cls)) * 100, 0 ), " %" ), collapse = " " ) )
      
      return( sum( as.integer(mapping) & as.integer( sim_list_entry ) ) )
    }
    adv <<- 0
    
    res_table = data.frame(
      
      "CL_name" = names(sim_list),
      "Intersect" = as.double( unlist( lapply( sim_list, FUN = match_fp ) ) ),
      "All_mutations" = unlist( lapply( sim_list, FUN = sum ) ),
      "Passed_treshold" = rep(F, nr_cls )
    )
    
    res_table = res_table[ order( res_table$Intersect, decreasing = T),  ]
    res_table$Passed_threshold[  res_table$Intersect >= 2] = T
    res_table$Passed_threshold[  res_table$Intersect < 2]  = F
    
    if ( dim(res_table[ res_table$Passed_threshold,])[1] >= 1 ){
      
      print( paste0( "Candidate(s): ", paste0( (res_table$CL_name[ res_table$Passed_threshold  ] ), collapse = "," ) )  )
      
      library("stringr")
      
      res_table[ res_table$Passed_threshold,]
      
      if ( output_path == "" )
        output_path = paste0( vcf_file_path, ".identification.tab" )
      
      print( paste0("Storing information in table: ",output_path ) )
      
    } else{
      print( paste0("No CL mutational fingerprint with sufficient similarity found." ) )
    }
    
    write.table( file = output_path_panel, res_table, sep ="\t", row.names = F, quote = F  )
  }
}}