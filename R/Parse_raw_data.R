### parse files
#' Parses data into r list variable
#' @param cosmic_genotype_file The path to the cosmic DNA genotype data file. Ensure that the right reference genome is used
#' @param ccle_genotype_file The path to the ccle DNA genotype data file. Ensure that the right reference genome is used
#' @param ref_gen Reference genome version
#' @param distinct_mode Should the mutational weights be calculated for all panels together or each for itelf? Recommendation: Seperately
#' @export
initiate_canonical_databases = function(
    cosmic_file = "CosmicCLP_MutantExport.tsv",
    ccle_file = "CCLE_hybrid_capture1650_hg19_NoCommonSNPs_CDS_2012.05.07.maf",
    ref_gen = "GRCH37",
    distinct_mode = TRUE
  ){
  
  require( "dplyr", quietly = TRUE, warn.conflicts = FALSE)
  require( "stringr", quietly = TRUE, warn.conflicts = FALSE)
  
  print( c( "Reference genome: ", ref_gen )  )
  
  ### pre processing
  
  package_path    = system.file("", package="Uniquorn")
  database_path   =  paste( package_path, "uniquorn_distinct_panels_db.sqlite3", sep ="/" )
  
  if (!distinct_mode)
    database_path   =  paste( package_path, "uniquorn_non_distinct_panels_db.sqlite3", sep ="/" )
  
  database_default_path =  paste( package_path, "uniquorn_db_default.sqlite3", sep ="/" )
  sim_list = as.data.frame( tbl( src_sqlite( database_default_path ), "sim_list_df" ), n = -1 )
  sim_list = sim_list[, which( colnames(sim_list) != "Ref_Gen"  ) ]
  sim_list = sim_list[, which( colnames(sim_list) != "Weight"  ) ]

  parse_files = c()
  
  if (file.exists(cosmic_genotype_file)){
      
      print( c( "Found CoSMIC: ", file.exists(cosmic_genotype_file) )  )
      sim_list = parse_cosmic_genotype_data( cosmic_genotype_file, sim_list )
      parse_files = c(parse_files, cosmic_genotype_file)
  }
  
  if (file.exists(ccle_genotype_file)){
      
      print( c( "Found CCLE: ", file.exists( ccle_genotype_file ) )  )
      sim_list = parse_ccle_genotype_data( ccle_genotype_file, sim_list )
      parse_files = c(parse_files, ccle_genotype_file)
  }
  
  if (length(parse_files) == 0)
      stop("Did not find CCLE & CoSMIC CLP file! Aborting.")
  
  # overwrite existing db
  if (file.exists(database_path))
    file.remove( database_path )

   print("Started pre-calculations")
  
  if ( exists("sim_list_stats"))
    rm( sim_list_stats )
  
  print("Finished parsing, aggregating over parsed Cancer Cell Line data")
  
  list_of_cls = unique( sim_list$CL )
  panels = sapply( list_of_cls, FUN = str_split, "_"  )
  panels = as.character(unique( as.character( sapply( panels, FUN = tail, 1) ) ))
  
  if (!distinct_mode){
    panels = paste0( c(panels), collapse ="|"  )
    database_path =  paste( package_path, "uniquorn_non_distinct_panels_db.sqlite3", sep ="/" )
  }
  
  print( paste( "Distinguishing between panels:",paste0( c(panels), collapse = ", "), sep = " ") )
  
  for (panel in panels) {
  
    print(panel)
    
    sim_list_panel   = sim_list[ grepl( panel, sim_list$CL) , ]
    member_var_panel = rep( 1, dim(sim_list_panel)[1] )
    
    sim_list_stats_panel = aggregate( member_var_panel , by = list( sim_list_panel$CL ), FUN = sum )
    colnames(sim_list_stats_panel) = c( "CL", "Count" )
    
    print("Aggregating over mutational frequency to obtain mutational weight")
        
    weights_panel = aggregate( member_var_panel , by = list( sim_list_panel$Fingerprint ), FUN = sum )
    weights_panel$x = 1.0 / as.double( weights_panel$x )
    
    mapping_panel = match( as.character( sim_list_panel$Fingerprint ), as.character( weights_panel$Group.1) )
    sim_list_panel = cbind( sim_list_panel, weights_panel$x[mapping_panel] )
    colnames( sim_list_panel )[3] = "Weight"
    
    # calculate weights
    
    aggregation_all_panel = stats::aggregate( 
      x  = as.double( sim_list_panel$Weight ),
      by = list( as.character( sim_list_panel$CL ) ),
      FUN = sum
    )
    
    mapping_agg_stats_panel = which( aggregation_all_panel$Group.1 %in% sim_list_stats_panel[,1], arr.ind = T  )
    sim_list_stats_panel = cbind( sim_list_stats_panel, aggregation_all_panel$x[mapping_agg_stats_panel] )
    
    #print("Finished aggregating, writing to database")
    
    Ref_Gen = rep(ref_gen, dim(sim_list_panel)[1]  )
    sim_list_panel = cbind( sim_list_panel, Ref_Gen )
    Ref_Gen = rep( ref_gen, dim(sim_list_stats_panel)[1]  )
    sim_list_stats_panel = cbind( sim_list_stats_panel, Ref_Gen )
    colnames( sim_list_stats_panel ) = c( "CL","Count","All_weights","Ref_Gen" )
    
    if(! exists("sim_list_global"))
      sim_list_global <<- sim_list[0,]
    
    sim_list_global = rbind(sim_list_global,sim_list_panel)
    
    if(! exists("sim_list_stats_global"))
      sim_list_stats_global <<- sim_list_stats_panel[0,]
    
    sim_list_stats_global = rbind( sim_list_stats_global, sim_list_stats_panel  )
  }
  
  print("Finished aggregating, saving to database")
  
  uni_db            = src_sqlite( path = database_path, create = T )
  sim_list_df       = tbl_df( sim_list_global )
  sim_list_stats_df = tbl_df( sim_list_stats_global )
  
  copy_to( uni_db, sim_list_df, temporary = F, 
    indexes = list(
      "Fingerprint",
      "CL",
      "Weight",
      "Ref_Gen"
    )
  )
  
  copy_to( uni_db, sim_list_stats_df, temporary = F,
    indexes = list(
      "CL",
      "Count",
      "All_weights",
      "Ref_Gen"
    )
  )
  
  print ("Initialization of Uniquorn DB finished")
}