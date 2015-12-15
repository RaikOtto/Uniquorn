#' Parses the vcf file and predicts the identity of the sample
#' @export
identify_vcf_file = function( 
  vcf_file_path, 
  output_path = "", 
  panels = c("CELLMINER","CCLE","COSMIC"), 
  type = "non_unique", 
  ref_gen = "grch37" ){
  # types = c("unique","non_unique")
  library("stringr")
  
  print( paste0( "Creating fingerprint from VCF file ", vcf_file_path  ) )
  vcf_fingerprint = parse_vcf_file( vcf_file_path )
  
  sim_list_store_mat = matrix( list(), ncol = length( type ), nrow = length( panels ) , dimnames = list( panels, type ) )
  sim_list_stats_store_mat = matrix( list(), ncol = length( type ), nrow = length( panels ) , dimnames = list( panels, type ) )
  
  if ( output_path == ""  ){
    
    output_file = paste( vcf_file_path, "ident.tab", sep ="_")
    
  } else {
    
    output_file = paste( paste0( output_path, tail( unlist(str_split(vcf_file_path, "/" ) ), 1 )  ), ".tab", sep ="_ident")
    
  }
  
  res_table <<- data.frame(
    
    "CL" = as.character(),
    "CL_source" = as.character(),
    "Found_muts_abs" = as.character(),
    "Count_mutations_abs" = as.character(),
    "Found_muts_rel" = as.character(),
    "Found_muts_weighted" = as.character(),
    "Found_muts_weighted_rel" = as.character(),
    "Count_mutations_weighted" = as.character(),
    "Passed_threshold" = as.character()
  )

  for( panel in panels ){
    
    if (type == "unique"){
        
      sim_list_file = paste( system.file("", package = "Uniquorn"), paste0( c( type,
      "_parsed_DB_mut_labels_",
      panel,".tab"),collapse= ""),sep = "/")

      sim_list_stats_file = paste( system.file("", package = "Uniquorn"),paste0( c( type, 
      "_parsed_DB_mut_labels_stats",
      panel,".tab"),collapse= ""),sep = "/")

    } else {
      
      sim_list_file = paste( system.file("", package = "Uniquorn"), paste0( c( type,
      "_parsed_DB_weighted_mut_labels_",
      panel,".tab"),collapse= ""),sep = "/")
      
      sim_list_stats_file = paste( system.file("", package = "Uniquorn"),paste0( c( type, 
      "_parsed_DB_weighted_mut_labels_stats_",
      panel,".tab"),collapse= ""),sep = "/")
    }

    print( paste0( "Loading similarity data from file ",  sim_list_file )  )
    
    if ( length( sim_list_store_mat[[ panel, type ]])  != 0 ){
        
      sim_list       = sim_list_store_mat[[ panel, type ]]
      sim_list_stats = sim_list_stats_store_mat[[ panel, type ]]
      
    } else {
        
        sim_list       = read.table( sim_list_file, sep = "\t", header = T)
        sim_list_stats = read.table( sim_list_stats_file, sep = "\t", header = T)
        
        sim_list_store_mat[[ panel, type ]] = sim_list
        sim_list_stats_store_mat[[ panel, type ]] = sim_list_stats
    }

    print( panel  )
    
    list_of_cls       = sort( unique( sim_list$CL ) )
    nr_cls            = length( list_of_cls  ) # amount cls
    found_mut_mapping = which( sim_list$Fingerprint %in% vcf_fingerprint ) # mapping
    cl_match_stats    = match( sim_list_stats$CL, list_of_cls  ) # mapping
    
    candidate_hits_abs = table( sim_list$CL[ found_mut_mapping ] )
    candidate_hits_abs = as.integer( candidate_hits_abs[  match( names( candidate_hits_abs) , list_of_cls  )  ] )
    candidate_hits_rel = round( candidate_hits_abs / sim_list_stats$Count[ cl_match_stats ], 5 ) * 100
    
    # weights
    cl_weight     = rep( 0.0, nr_cls )
    cl_weight_rel = rep( 0.0, nr_cls )
    all_weighted  = rep( 0.0, nr_cls )
    
    # threshold non weighted
    passed_threshold_vec = rep( F, nr_cls )
    passed_threshold_vec[ ( candidate_hits_abs >= 3 ) & ( candidate_hits_rel >= 2 ) ] = T
    
    if (type != "unique"){
      
       # aggregate over weights & CL identifier
        
        aggregation = stats::aggregate(
          x  = as.double( sim_list$Weight[ found_mut_mapping  ] ),
          by = list( as.character( sim_list$CL[ found_mut_mapping ] )  ),
          FUN = sum
        )
        
        aggregation_all = stats::aggregate( 
          x  = as.double( sim_list$Weight ),
          by = list( as.character( sim_list$CL ) ),
          FUN = sum
        )
        
        aggregation_all = aggregation_all[ match( list_of_cls, aggregation_all[, 1]  )  ,]
        all_weighted    = aggregation_all[ match( list_of_cls, aggregation_all[, 1]  ), 2]
        aggregation_match = match( aggregation[,1], as.character( list_of_cls )  )
        
        cl_weight[ aggregation_match ] = aggregation[ , 2 ]
        cl_weight_rel = round( as.double( cl_weight ) / as.double( aggregation_all[ ,2 ] ) , 3 ) * 100
        
        passed_threshold_vec = rep( F, nr_cls )
        passed_threshold_vec[ ( cl_weight >= 10 ) & ( candidate_hits_abs >= 2 ) ] = T

    }
    
    res_table = data.frame( 
      "CL"                      = c( as.character( res_table$CL ) , as.character( list_of_cls ) ),
      "CL_source"               = c( as.character( res_table$CL_source), as.character( rep(panel, nr_cls) ) ),
      "Found_muts_abs"          = c( as.character( res_table$Found_muts_abs), as.character( candidate_hits_abs ) ),
      "Count_mutations_abs"     = c( as.character( res_table$Count_mutations_abs), as.character(  sim_list_stats$Count[cl_match_stats] ) ),
      "Found_muts_rel"          = c( as.character( res_table$Found_muts_rel), as.character(  candidate_hits_rel ) ),
      "Found_muts_weighted"     = c( as.character( res_table$Found_muts_weighted ),cl_weight ),
      "Count_mutations_weighted"= c( as.character( res_table$Count_mutations_weighted), as.character(  all_weighted ) ),
      "Found_muts_weighted_rel" = c( as.character( res_table$Found_muts_weighted_rel ),cl_weight_rel ),
      "Passed_threshold"        = c( as.character( res_table$Passed_threshold), as.character( passed_threshold_vec ) )
    )

    # don't sort me bro!
    
    
  }
  
  if (type != "unique"){
    res_table = res_table[ order( as.double( res_table$Found_muts_weighted_rel ), decreasing = T),  ]
  } else {
    res_table = res_table[ order( as.double( res_table$Found_muts_rel ), decreasing = T),  ]
  }
  
  if ( dim( res_table[ res_table$Passed_threshold ,])[1] == 1 ){
    
    print( paste0( "Candidate(s): ", paste0( ( res_table$CL_name[ res_table$Passed_treshold  ] ), collapse = "," ) )  )
    
  } else if ( dim( res_table[ res_table$Passed_threshold ,])[1] > 1 ){
    
    print( paste0( "Candidates: ", paste0( ( res_table$CL_name[ res_table$Passed_treshold  ] )[1:2], collapse = "," ) )  )
    
  } else {
    
    print( paste0("No CL mutational fingerprint with sufficient similarity found." ) )
  }  
  
  print( paste0("Storing information in table: ",output_file ) )
  write.table( res_table, output_file, sep ="\t", row.names = F, quote = F  )
}