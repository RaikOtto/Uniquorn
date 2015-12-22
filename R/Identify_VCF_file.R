#' Parses the vcf file and predicts the identity of the sample
#' @export
identify_vcf_file = function( 
  vcf_file_path,
  output_path = "",
  ref_gen = "hg19" ){
  library( "stringr" )
  
  message( paste0("Assuming reference genome ", ref_gen) )
  
  ### pre processing
  
  path_to_python  = paste( system.file("", package="Uniquorn"),"pre_compute_raw_data.py", sep ="/")
  
  uni_db_path =  paste( db_folder, paste0( c( ref_gen , "uniquorn_db.sqlite3"), collapse = "_"), sep ="/" )
  
  # reading file
  print( paste0( "Creating fingerprint from VCF file ", vcf_file_path  ) )
  vcf_fingerprint = parse_vcf_file( vcf_file_path )
  
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
    "Passed_threshold" = as.character(),
    "Found_muts_weighted" = as.character(),
    "Found_muts_weighted_rel" = as.character(),
    "Count_mutations_weighted" = as.character(),
    "Passed_threshold_weighted" = as.character()
  )
    
    if ( ! file.exists( uni_db_path ) ){
      
      message( paste0( "Did not find database for reference genome : ", ref_gen ) )
      
    } else {

      print( paste0( "Loading similarity database for reference genome ",  ref_gen )  )
      
      
      
      list_of_cls       = sort( unique( sim_list$CL ) )
      nr_cls            = length( list_of_cls  ) # amount cls
      found_mut_mapping = which( sim_list$Fingerprint %in% vcf_fingerprint ) # mapping
      cl_match_stats    = match( list_of_cls, sim_list_stats$CL   ) # mapping
      
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
      passed_threshold_weighted = rep( "", nr_cls )
      
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
      
      passed_threshold_weighted = rep( F, nr_cls )
      passed_threshold_weighted[ cl_weight_rel >= 10.0 ] = T
      
      ouput_cl_names = str_replace( as.character( list_of_cls ), pattern = paste0( "_", panel  ), replacement = "" )
      
      res_table = data.frame( 
        "CL"                       = c( as.character( res_table$CL ) , ouput_cl_names ),
        "CL_source"                = c( as.character( res_table$CL_source), as.character( rep(panel, nr_cls) ) ),
        "Found_muts_abs"           = c( as.character( res_table$Found_muts_abs), as.character( candidate_hits_abs ) ),
        "Count_mutations_abs"      = c( as.character( res_table$Count_mutations_abs), as.character(  sim_list_stats$Count[cl_match_stats] ) ),
        "Found_muts_rel"           = c( as.character( res_table$Found_muts_rel), as.character(  candidate_hits_rel ) ),
        "Passed_threshold"         = c( as.character( res_table$Passed_threshold), as.character( passed_threshold_vec ) ),
        "Found_muts_weighted"      = c( as.character( res_table$Found_muts_weighted ),cl_weight ),
        "Count_mutations_weighted" = c( as.character( res_table$Count_mutations_weighted), as.character(  all_weighted ) ),
        "Found_muts_weighted_rel"  = c( as.character( res_table$Found_muts_weighted_rel ),cl_weight_rel ),
        "Passed_threshold_weighted"= c( as.character( res_table$Passed_threshold_weighted), as.character( passed_threshold_weighted ) )
      )
      }
  
  res_table = res_table[ order( as.double( res_table$Found_muts_weighted_rel ), decreasing = T),  ]
  
  print( paste0( "Candidate(s): ", paste0( ( unique( as.character( res_table$CL )[ res_table$Passed_threshold_weighted == T  ]) ), collapse = "," ) )  )
  
  print( paste0("Storing information in table: ",output_file ) )
  write.table( res_table, output_file, sep ="\t", row.names = F, quote = F  )
}