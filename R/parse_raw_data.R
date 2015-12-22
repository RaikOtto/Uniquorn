### parse files
#' Parses data into r list variable
#' @export
initiate_canonical_databases = function(
    cosmic_genotype_file = "CosmicCLP_MutantExport.tsv",
    cellminer_genotype_file = 'DNA__Exome_Seq_none.txt',
    ccle_genotype_file = "CCLE_hybrid_capture1650_hg19_allVariants_2012.05.07.maf",
    overwrite_old_db = F,
    ref_gen = "hg19"
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
  
  uni_db_path       =  paste( db_folder, paste0( c( ref_gen , "uniquorn_db.sqlite3"), collapse = "_"), sep ="/" )
  uni_db_stats_path =  paste( db_folder, paste0( c( ref_gen , "uniquorn_db_stats.sqlite3"), collapse = "_"), sep ="/" )
  
  if ( overwrite_old_db){
    if ( file.exists(uni_db_path))
      file.remove(uni_db_path)
    if ( file.exists(uni_db_stats_path))
      file.remove(uni_db_stats_path)
  }
  
  uni_db       = src_sqlite( uni_db_path, create = T)
  uni_db_stats = src_sqlite( uni_db_stats_path, create = T)
  
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
  
  message("Parsing data finished, loading sqlite database")
  
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
  
  uni_db       = src_sqlite( uni_db_path, create = T )
  
  sim_list_df     = tbl_df(sim_list)
  sim_list_sqlite = copy_to( uni_db, sim_list_df, temporary = F, indexes = list(
    "Fingerprint",
    "CL",
    "Weight"
    )
  )
  
  uni_db_stats = src_sqlite( uni_db_stats_path, create = T)
  
  sim_list_stats_df     = tbl_df(sim_list_stats)
  sim_list_stats_sqlite = copy_to( 
    uni_db_stats, 
    sim_list_stats_df, 
    temporary = FALSE,
    indexes = list(
      "CL",
      "Count"
    )
  )
  
  print ("Finished initializing the canonical Cancer Cell line trainingssets")
}