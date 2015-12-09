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
    
    fingerprint_stats = read.table( fingerprint_names_file, sep ="\t", header= F )
    all_fingerprints_in_cl_set = fingerprint_stats[,1]
    all_weights_muts_in_cl_set = fingerprint_stats[,2]
    
    print( "Finished loading similarity data. Mapping vcf's fingerprint to all contained fingerprints"  )
    
    mapping = match( all_fingerprints_in_cl_set, vcf_fingerprint, nomatch = 0 )
    
    adv = 0
    nr_cls = length( sim_list )
    
    match_fp = function( sim_list_entry ){
      
      stat = round( (adv / as.double(nr_cls)) * 100, 0 )
      adv <<- adv + 1
      
      if ( stat != round( (adv / as.double(nr_cls)) * 100, 0 ) )
        print( paste0( c( panel,  round( (adv / as.double(nr_cls)) * 100, 0 ), " %" ), collapse = " " ) )
      
      res_map = ( as.integer(mapping) & as.integer( sim_list_entry ) ) 
      res_score = sum( all_weights_muts_in_cl_set[ res_map ] )
      
      res_scores = c( sum(res_map), res_score )
      
      return( res_scores )
    }
    adv <<- 0
    
    res_intersect = as.data.frame( lapply( sim_list, FUN = match_fp ) )
    
    res_table = data.frame(
      
      "CL_name" = as.character( names(sim_list) ),
      "Intersect" = as.double( res_intersect[1,] ),
      "Intersect_weighted" = as.double( res_intersect[2,] ),
      "All_mutations" = as.integer( unlist( lapply( sim_list, FUN = sum ) ) ),
      "Passed_treshold" = rep( F, nr_cls ),
      "Passed_treshold_weighted" = rep( F, nr_cls )
    )
    
    res_table = res_table[ order( as.double( res_table$Intersect_weighted ), decreasing = T),  ]
    
    res_table$Passed_treshold[  
      ( as.double(res_table$Intersect) >= 2.0) & ( ( as.double(res_table$Intersect) / as.double(res_table$All_mutations)) >= .01  )
    ] = T
    res_table$Passed_treshold[ is.na(res_table$Passed_treshold) ]  = F
    
    res_table$Passed_treshold_weighted[  
      ( as.double(res_table$Intersect_weighted) >= 10.0) & ( ( as.double(res_table$Intersect_weighted) / as.double(res_table$All_mutations)) >= .03  )
      ] = T
    res_table$Passed_treshold_weighted[ is.na(res_table$Passed_treshold_weighted) ]  = F
    
    
    if ( dim( res_table[ res_table$Passed_treshold_weighted ,])[1] >= 1 ){
      
      print( paste0( "Candidate(s): ", paste0( (res_table$CL_name[ res_table$Passed_treshold_weighted  ] ), collapse = "," ) )  )
      
    } else{
      print( paste0("No CL mutational fingerprint with sufficient similarity found." ) )
    }
    
    print( paste0("Storing information in table: ",output_path_panel ) )
    write.table( file = output_path_panel, res_table, sep ="\t", row.names = F, quote = F  )
  }
}}