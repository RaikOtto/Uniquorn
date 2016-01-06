### parse files
#' Parses data into r list variable
#' @export
initiate_canonical_databases = function(
    cosmic_genotype_file = "CosmicCLP_MutantExport.tsv",
    ccle_genotype_file = "CCLE_hybrid_capture1650_hg19_allVariants_2012.05.07.maf",
    ref_gen = "HG19"
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
  db_folder       = system.file("", package="Uniquorn")
  
  uni_db_path       =  paste( db_folder, "uniquorn_db.sqlite3", sep ="/" )
  uni_db_default_path =  paste( db_folder, "uniquorn_db_default.sqlite3", sep ="/" )

  if ( file.exists(uni_db_path) )
    file.copy( uni_db_default_path, uni_db_path)

  # python parser
  
  print("Started parsing")
  
  command_line = str_c( 
    c(  
      'python',     path_to_python,
      "-ccle ",     ccle_genotype_file,
      "-cosmic ",   cosmic_genotype_file,
      "-cellminer", cellminer_genotype_file,
      "-o_db_path", db_folder
    ),
    collapse = " "
  )
  
  system( command_line, ignore.stdout = F, intern = F )
  
  if ( exists("sim_list"))
    rm( sim_list )
  if ( exists("sim_list_stats"))
    rm( sim_list_stats )
  
  for( parse_file in parse_files ){
   
    if (parse_file == cosmic_genotype_file){
      
      panel = "cosmic"
    
    } else if (parse_file == ccle_genotype_file){
      
      panel = "ccle"
      
    }
    
    print( paste( "Parsing: ", panel ), sep =" "  )
    
    sim_list_file = paste0( c( db_folder, "/", "Fingerprint_",       panel, ".tab" ), collapse = "" )
    
    if (! exists("sim_list")){
      sim_list = read.table( sim_list_file, sep = "\t", header = T)
    } else {
      sim_list = rbind(sim_list, read.table( sim_list_file, sep = "\t", header = T))
    }
    
    #file.remove(sim_list_file)
    #file.remove(sim_list_stats_file)
  }
  
  print("Finished parsing, aggregating over parsed Cancer Cell Line data")
  member_var = rep( 1, dim(sim_list)[1] )
  sim_list_stats = aggregate( member_var , by = list( sim_list$CL ), FUN = sum )
  colnames(sim_list_stats) = c( "CL", "Count" )
  
  print("Aggregating over mutational frequency to obtain mutational weight")

  weights = aggregate( member_var , by = list( sim_list$Fingerprint ), FUN = sum )
  weights$x = 1.0 / as.double( weights$x )
  
  mapping = match( as.character( sim_list$Fingerprint ), as.character( weights$Group.1) )
  sim_list = cbind( sim_list, weights$x[mapping] )
  colnames( sim_list )[3] = "Weight"
  
  print("Finished aggregating, writing to database")
  
  Ref_Gen = rep(ref_gen, dim(sim_list)[1]  )
  sim_list = cbind( sim_list, Ref_Gen )
  Ref_Gen = rep( ref_gen, dim(sim_list_stats)[1]  )
  sim_list_stats = cbind( sim_list_stats, Ref_Gen )
  
  uni_db            = src_sqlite( uni_db_path, create = T )
  sim_list_df       = tbl_df( sim_list )
  sim_list_stats_df = tbl_df( sim_list_stats )
  
  copy_to( uni_db, sim_list_df, temporary = F, indexes = list(
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
      "Ref_Gen"
    )
  )
  
  print ("Initialization of Uniquorn DB finished")
}