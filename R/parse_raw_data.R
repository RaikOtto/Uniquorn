### parse files
#' Parses data into r list variable
#' @export
initiate_canonical_databases = function(
    cosmic_genotype_file = "CosmicCLP_MutantExport.tsv",
    cellminer_genotype_file = 'DNA__Exome_Seq_none.txt',
    ccle_genotype_file = "CCLE_hybrid_capture1650_hg19_allVariants_2012.05.07.maf",
    ref_gen = "HG19"
  ){
  
  suppressPackageStartupMessages(library("plyr"))
  suppressPackageStartupMessages(library("dplyr"))
  suppressPackageStartupMessages(library("stringr"))
  
  parse_files = c()
  
  if (file.exists(cosmic_genotype_file)){
    
    print( c( "Found CoSMIC: ", file.exists(ccle_genotype_file) )  )
    parse_files = c(parse_files, cosmic_genotype_file)
  }
  
  if (file.exists(ccle_genotype_file)){
    
    print( c( "Found CCLE: ", file.exists(ccle_genotype_file) )  )
    parse_files = c(parse_files, ccle_genotype_file)
  }
  
  if (file.exists(cellminer_genotype_file)){
    
    print( c( "Found CellMiner: ", file.exists(cellminer_genotype_file) )  )
    parse_files = c(parse_files, cellminer_genotype_file)
  }
  
  print( c( "Reference genome: ", ref_gen )  )
  
  ### pre processing
  
  path_to_python  = paste( system.file("", package="Uniquorn"),"pre_compute_raw_data.py", sep ="/")
  db_folder       = system.file("", package="Uniquorn")
  
  uni_db_path       =  paste( db_folder, "uniquorn_db.sqlite3", sep ="/" )

  if ( file.exists(uni_db_path) )
    file.remove(uni_db_path)

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
  print("Finished parsing of raw data, transforming data")
  
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
      
    } else if (parse_file == cellminer_genotype_file){
      
      panel = "cellminer"
    }
    
    sim_list_file       = paste0( c( db_folder, "/", "Fingerprint_",       panel, ".tab" ), collapse = "" )
    sim_list_stats_file = paste0( c( db_folder, "/", "Fingerprint_stats_", panel, ".tab" ), collapse = "" )
    
    if (! exists("sim_list")){
      sim_list = read.table( sim_list_file, sep = "\t", header = T)
    } else {
      sim_list = rbind(sim_list, read.table( sim_list_file, sep = "\t", header = T))
    }
    
    if (! exists("sim_list_stats")){
      sim_list_stats = read.table( sim_list_stats_file, sep = "\t", header = T)
    } else {
      sim_list_stats = rbind(sim_list_stats, read.table( sim_list_stats_file, sep = "\t", header = T))
    }
    
    file.remove(sim_list_file)
    file.remove(sim_list_stats_file)
  }
  
  Ref_Gen = rep(ref_gen, dim(sim_list)[1]  )
  sim_list = cbind( sim_list, Ref_Gen )
  Ref_Gen = rep( ref_gen, dim(sim_list_stats)[1]  )
  sim_list_stats = cbind( sim_list_stats, Ref_Gen )
  
  print("Finished transforming data, writing to database")
  
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
  
  print ("Finished initializing the canonical Cancer Cell line trainingssets")
}