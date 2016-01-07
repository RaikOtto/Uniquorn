#' Parses the vcf file and predicts the identity of the sample
#' @export
identify_vcf_file = function( 
  vcf_file,
  output_file = "",
  ref_gen = "HG19",
  similarity_threshold = 15.0,
  unique_mode = F){
  
  suppressPackageStartupMessages( library( "stringr" ) )
  suppressPackageStartupMessages( library( "dplyr" ) )
  suppressPackageStartupMessages( library( "plyr" ) )
  
  message( paste0("Assuming reference genome ", ref_gen) )
  
  ### pre processing
  
  package_path    = system.file("", package="Uniquorn")
  path_to_python  = paste( package_path,"pre_compute_raw_data.py", sep ="/")
  database_path     =  paste( package_path, "uniquorn_db.sqlite3", sep ="/" )
  
  # reading file
  vcf_fingerprint = parse_vcf_file( vcf_file )
  
  if ( output_file == ""  )
    output_file = paste( vcf_file, "uniquorn_ident.tab", sep ="_")
    
  if( ! file.exists( database_path ) ){
    
    database_path = paste( package_path, "uniquorn_db_default.sqlite3", sep ="/" )
    message("CCLE & CoSMIC CLP cancer cell line fingerprint NOT found, defaulting to 60 CellMiner cancer cell lines! 
            We strongly advise to add CCLE & CoSMIC, see readme.")
  }
    
  print( "Finished reading the VCF file, loading database" )
  
  sim_list = as.data.frame( tbl( src_sqlite( database_path ), "sim_list_df" ), n = -1 )
  
  sim_list = sim_list[ sim_list$Ref_Gen == ref_gen  ,]
  print(paste0( c("Found ", as.character( length( unique(sim_list$CL) ) ), " many CLs for reference genome ", ref_gen ), collapse = "" ) )
  
  sim_list_stats = as.data.frame( tbl( src_sqlite( database_path ), "sim_list_stats_df" ), n = -1 )
  sim_list_stats = sim_list_stats[ sim_list_stats$Ref_Gen == ref_gen  ,]
  
  print("Finished reading database, identifying CL")
  
  if ( unique_mode  ){
    
    print("Unique mode, only using mutations that are unique to cancer cell lines")
    
    sim_list = sim_list[ sim_list$Weight == 1 ,  ]
    sim_list_stats = aggregate( as.double( sim_list$Weight ), by = list( sim_list$CL), FUN = sum  )
    colnames( sim_list_stats ) = c("CL","Count")
  }
  
  list_of_cls       = unique( sim_list$CL )
  nr_cls            = length( list_of_cls  ) # amount cls
  
  found_mut_mapping = which( sim_list$Fingerprint %in% vcf_fingerprint ) # mapping

  ### unweighted scores
  
  candidate_hits_abs_all = rep(0, nr_cls)
  names(candidate_hits_abs_all) = list_of_cls
  
  candidate_hits_abs = aggregate( 
    rep(1, length(found_mut_mapping)),
    by = list(sim_list$CL[ found_mut_mapping ]), 
    FUN = sum
  )
  candidate_hits_abs_all[ which( list_of_cls %in% candidate_hits_abs$Group.1) ] = as.integer( candidate_hits_abs$x[ match( list_of_cls, candidate_hits_abs$Group.1, nomatch = 0 ) ] )
  
  cl_match_stats    = match( list_of_cls, sim_list_stats$CL, nomatch = 0 ) # mapping
  candidate_hits_rel = round( candidate_hits_abs_all / sim_list_stats$Count[ cl_match_stats ] * 100, 1 ) 
  
  ### weighted scores
  
  # weights
  
  cl_weight     = rep( 0.0, nr_cls )
  cl_weight_rel = rep( 0.0, nr_cls )
  all_weighted  = rep( 0.0, nr_cls )
  
  # threshold non weighted
  passed_threshold_vec = rep( F, nr_cls )
  passed_threshold_vec[ ( candidate_hits_abs_all >= 3 ) & ( candidate_hits_rel >= 2 ) ] = T
  passed_threshold_weighted = rep( "", nr_cls )
  
  # aggregate over weights & CL identifier
  
  aggregation = stats::aggregate(
    x  = as.double( sim_list$Weight[ found_mut_mapping  ] ),
    by = list( as.character( sim_list$CL[ found_mut_mapping ] )  ),
    FUN = sum
  )
  
  weight_all = sim_list_stats$All_weights[ match( aggregation$Group.1, sim_list_stats$CL )  ]
  cl_weight_rel = round( as.double( aggregation$x ) / as.double( weight_all ) , 3 ) * 100
  mapping_to_cls = match( list_of_cls, aggregation$Group.1 , nomatch = 0  )
  
  res_cl_weighted = rep(0, nr_cls)
  names(res_cl_weighted) = list_of_cls
  res_cl_weighted[ names(res_cl_weighted) %in% aggregation$Group.1  ] = aggregation$x[ mapping_to_cls ]
  stats_all_weight = sim_list_stats$All_weights[ match( list_of_cls, sim_list_stats$CL  ) ]
  
  res_res_cl_weighted = round( as.double(res_cl_weighted  ) / stats_all_weight * 100, 1 )
  res_cl_weighted = round(res_cl_weighted, 0)
  
  # treshold
  
  passed_threshold_weighted = rep( F, nr_cls )
  passed_threshold_weighted[ (res_res_cl_weighted >= 15.0) ] = TRUE
  
  if( unique_mode ){
    
    passed_threshold_weighted = rep( F, nr_cls )
    passed_threshold_weighted[ ( candidate_hits_abs_all >= 3.0) & ( candidate_hits_rel >= .05 ) ] = T
  }

  output_cl_names = str_replace( list_of_cls, pattern = "_CCLE|_COSMIC|_CELLMINER", replacement = "" )
  panel_vec = rep("", length( output_cl_names ))
  panel_vec[ str_detect( list_of_cls, "_CCLE" ) ] = "CCLE"
  panel_vec[ str_detect( list_of_cls, "_COSMIC" ) ] = "COSMIC"
  panel_vec[ str_detect( list_of_cls, "_CELLMINER" ) ] = "CELLMINER"
  
  res_table = data.frame( 
    "CL"                       = output_cl_names,
    "CL_source"                = panel_vec,
    "Found_muts_abs"           = as.character( candidate_hits_abs_all ),
    "Count_mutations_abs"      = as.character(  sim_list_stats$Count[ cl_match_stats ] ),
    "Found_muts_rel"           = as.character(  candidate_hits_rel ),
    "Found_muts_weighted"      = as.character( res_cl_weighted ),
    "Count_mutations_weighted" = as.character( round( stats_all_weight, 0 ) ),
    "Found_muts_weighted_rel"  = as.character( res_res_cl_weighted ),
    "Passed_threshold"         = as.character( passed_threshold_weighted )
  )
  
  res_table = res_table[ order( as.double( as.character( res_table$Found_muts_weighted_rel) ), decreasing = T),  ]
  
  print( paste0( "Candidate(s): ", paste0( ( unique( as.character( res_table$CL )[ res_table$Passed_threshold == T  ]) ), collapse = "," ) )  )
  
  print( paste0("Storing information in table: ",output_file ) )
  write.table( res_table, output_file, sep ="\t", row.names = F, quote = F  )
}