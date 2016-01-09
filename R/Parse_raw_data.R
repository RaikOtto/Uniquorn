### parse files
#' Parses data into r list variable
#' @export
initiate_canonical_databases = function(
    cosmic_genotype_file = "CosmicCLP_MutantExport.tsv",
    ccle_genotype_file = "CCLE_hybrid_capture1650_hg19_allVariants_2012.05.07.maf",
    ref_gen = "GRCH37",
    distinct_mode = TRUE
  ){
  
  suppressPackageStartupMessages(library("plyr"))
  suppressPackageStartupMessages(library("dplyr"))
  suppressPackageStartupMessages(library("stringr"))
  
  parse_files = c()
  
  if (file.exists(cosmic_genotype_file)){
    
    print( c( "Found CoSMIC: ", file.exists(cosmic_genotype_file) )  )
    parse_files = c(parse_files, cosmic_genotype_file)
  }
  
  if (file.exists(ccle_genotype_file)){
    
    print( c( "Found CCLE: ", file.exists( ccle_genotype_file ) )  )
    parse_files = c(parse_files, ccle_genotype_file)
  }

  print( c( "Reference genome: ", ref_gen )  )
  
  ### pre processing
  
  path_to_python  = paste( system.file("", package="Uniquorn"),"pre_compute_raw_data.py", sep ="/")
  package_path    = system.file("", package="Uniquorn")
  
  database_path   =  paste( package_path, "uniquorn_distinct_panels_db.sqlite3", sep ="/" )
  
  if (!distinct_mode)
    database_path   =  paste( package_path, "uniquorn_non_distinct_panels_db.sqlite3", sep ="/" )
  
  database_default_path =  paste( package_path, "uniquorn_db_default.sqlite3", sep ="/" )
  sim_list_default = as.data.frame( tbl( src_sqlite( database_default_path ), "sim_list_df" ), n = -1 )
  sim_list_default = sim_list_default[, which( colnames(sim_list_default) != "Ref_Gen"  ) ]
  sim_list_default = sim_list_default[, which( colnames(sim_list_default) != "Weight"  ) ]

  # overwrite existing db
  if (file.exists(database_path))
    file.remove( database_path )

  # python parser

  print("Started pre-calculations")

  command_line = str_c( 
    c(  
      'python',     path_to_python,
      "-ccle ",     ccle_genotype_file,
      "-cosmic ",   cosmic_genotype_file,
      "-o_db_path", package_path
    ),
    collapse = " "
  )
  
  system( command_line, ignore.stdout = F, intern = F )
  
  if ( exists("sim_list_stats"))
    rm( sim_list_stats )
  
  for( parse_file in parse_files ){
   
    if (parse_file == cosmic_genotype_file){
      
      panel = "cosmic"
    
    } else if (parse_file == ccle_genotype_file){
      
      panel = "ccle"
      
    }
    
    print( paste( "Parsing: ", panel ), sep =" "  )
    
    sim_list_file = paste0( c( package_path, "/", "Fingerprint_",       panel, ".tab" ), collapse = "" )
    sim_list_default      = rbind( sim_list_default, read.table( sim_list_file, sep = "\t", header = T))
    
    #file.remove(sim_list_file)
    #file.remove(sim_list_stats_file)
  }
  
  list_of_cls = unique( sim_list_default$CL )
  panels = sapply( list_of_cls, FUN = str_split, "_"  )
  panels = as.character(unique( as.character( sapply( panels, FUN = tail, 1) ) ))
  
  print("Finished parsing, aggregating over parsed Cancer Cell Line data")
  print( paste( "Distinguishing between panels:",paste0( c(panels), collapse = ", "), sep = " ") )
  
  if (!distinct_mode){
    panels = paste0( c(panels), collapse ="|"  )
    database_path =  paste( package_path, "uniquorn_non_distinct_panels_db.sqlite3", sep ="/" )
  }
  
  for (panel in panels) {
  
    print(panel)
    
    sim_list_panel = sim_list_default[ grepl( panel, sim_list_default$CL) , ]
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
      sim_list_global <<- sim_list_default[0,]
    
    sim_list_global = rbind(sim_list_global,sim_list_panel)
    
    if(! exists("sim_list_stats_global")){
      sim_list_stats_global <<- sim_list_stats_panel[0,]

      sim_list_stats_global = rbind( sim_list_stats_global, sim_list_stats_panel  )
    }
  }
  
  uni_db            = src_sqlite( database_path, create = T )
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