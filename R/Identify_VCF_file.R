#' Parses the vcf file and predicts the identity of the sample
#' @export
identify_vcf_file = function( 
  vcf_file,
  output_file = "",
  ref_gen = "HG19",
  unique_mode = F){
  suppressPackageStartupMessages( library( "stringr" ) )
  suppressPackageStartupMessages( library( "dplyr" ) )
  suppressPackageStartupMessages( library( "plyr" ) )
  
  message( paste0("Assuming reference genome ", ref_gen) )
  
  ### pre processing
  
  package_path    = system.file("", package="Uniquorn")
  path_to_python  = paste( package_path,"pre_compute_raw_data.py", sep ="/")
  uni_db_path     =  paste( package_path, "uniquorn_db.sqlite3", sep ="/" )
  
  # reading file
  vcf_fingerprint = parse_vcf_file( vcf_file )
  
  if ( output_file == ""  )
    output_file = paste( vcf_file, "uniquorn_ident.tab", sep ="_")
    
  if( ! file.exists( uni_db_path ) ){
    
    uni_db_path = paste( package_path, "uniquorn_db_default.sqlite3", sep ="/" )
    message("CCLE & CoSMIC CLP cancer cell line fingerprint NOT found, defaulting to 65 CellMiner cancer cell lines! We strongly advise to add CCLE & CoSMIC, see readme.")
  }
    
  message( "Finished reading the VCF file, loading database with trainingsets" )
  
  sim_list = as.data.frame( tbl( src_sqlite( uni_db_path ), "sim_list_df" ), n = -1 )
  
  sim_list = sim_list[ sim_list$Ref_Gen == ref_gen  ,]
  print(paste0( c("Found ", as.character( length( unique(sim_list$CL) ) ), " many CLs for reference genome ", ref_gen ), collapse = "" ) )
  
  sim_list_stats = as.data.frame( tbl( src_sqlite( uni_db_path ), "sim_list_stats_df" ), n = -1 )
  sim_list_stats = sim_list_stats[ sim_list_stats$Ref_Gen == ref_gen  ,]
  
  if ( unique_mode  ){
    
    print("Unique mode, only using mutations that are unique to cancer cell lines")
    
    sim_list = sim_list[ sim_list$Weight == 1 ,  ]
    sim_list_stats = aggregate( as.double( sim_list$Weight ), by = list( sim_list$CL), FUN = sum  )
    colnames( sim_list_stats ) = c("CL","Count")
  }
  
  list_of_cls       = unique( sim_list$CL )
  nr_cls            = length( list_of_cls  ) # amount cls
  
  found_mut_mapping = match( vcf_fingerprint, sim_list$Fingerprint, nomatch = 0 ) # mapping
  found_mut_mapping = found_mut_mapping[ found_mut_mapping != 0 ]

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
  candidate_hits_rel = round( candidate_hits_abs_all / sim_list_stats$Count[ cl_match_stats ], 5 ) * 100
  
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
  
  aggregation_all = stats::aggregate( 
    x  = as.double( sim_list$Weight ),
    by = list( as.character( sim_list$CL ) ),
    FUN = sum
  )
  
  aggregation_all = aggregation_all[ match( list_of_cls, aggregation_all[, 1]  )  ,]
  all_weighted    = aggregation_all[ match( list_of_cls, aggregation_all[, 1]  ), 2]
  aggregation_match = match( aggregation[,1], as.character( list_of_cls )  )
  
  cl_weight[ aggregation_match ] = round(aggregation[ , 2 ],0)
  cl_weight_rel = round( as.double( cl_weight ) / as.double( aggregation_all[ ,2 ] ) , 3 ) * 100
  
  # match value
  
  match_value = round( as.double( cl_weight ) / as.double( candidate_hits_abs_all ), 2 ) 
  match_value[ is.na( match_value) ] = 0
  
  # treshold
  
  passed_threshold_weighted = rep( F, nr_cls )
  passed_threshold_weighted[ (cl_weight_rel >= 5.0) & ( match_value >= .3 ) ] = T
  
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
    "Match_value"              = as.character( match_value ),
    "Found_muts_weighted"      = as.character( cl_weight ),
    "Count_mutations_weighted" = as.character( round( all_weighted, 0 ) ),
    "Found_muts_weighted_rel"  = as.character( cl_weight_rel ),
    "Passed_threshold"         = as.character( passed_threshold_weighted )
  )
  
  res_table = res_table[ order( as.double( as.character( res_table$Found_muts_weighted_rel) ), decreasing = T),  ]
  
  print( paste0( "Candidate(s): ", paste0( ( unique( as.character( res_table$CL )[ res_table$Passed_threshold == T  ]) ), collapse = "," ) )  )
  
  print( paste0("Storing information in table: ",output_file ) )
  write.table( res_table, output_file, sep ="\t", row.names = F, quote = F  )
}