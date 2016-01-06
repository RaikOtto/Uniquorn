
#' Show all cancer cell line identifier present in the database for a selected reference Genome
#' @export
show_contained_cls = function( ref_gen = "HG19" ){

  suppressPackageStartupMessages(library("plyr"))
  suppressPackageStartupMessages(library("dplyr"))
  suppressPackageStartupMessages(library("stringr"))
    
  db_folder      = system.file("", package="Uniquorn")
  database_path  =  paste( db_folder, "uniquorn_db.sqlite3", sep ="/" )
  sim_list_stats = as.data.frame( tbl( src_sqlite( database_path ), "sim_list_stats_df" ), n = -1 )

  sim_list_stats = sim_list_stats[ sim_list_stats$Ref_Gen == ref_gen,  ]
  print( paste0( c("Found ", dim(sim_list_stats)[1], " many cancer cell lines fingerprints for reference genome ", ref_gen ), collapse = ""  )  )
  
  print( summary( sim_list_stats ) )
  
  return( sim_list_stats )  
}

#' Show all mutations present in the database for a selected reference Genome
#' @export
show_contained_mutations = function( ref_gen = "HG19" ){
  
  suppressPackageStartupMessages(library("plyr"))
  suppressPackageStartupMessages(library("dplyr"))
  suppressPackageStartupMessages(library("stringr"))
  
  db_folder      = system.file("", package="Uniquorn")
  database_path  =  paste( db_folder, "uniquorn_db.sqlite3", sep ="/" )
  sim_list = as.data.frame( tbl( src_sqlite( database_path ), "sim_list_df" ), n = -1 )
  
  sim_list       = sim_list[ sim_list$Ref_Gen == ref_gen,  ]
  print( paste0( c("Found ", dim(sim_list)[1], " many cancer cell lines associated mutations for reference genome ", ref_gen ), collapse = ""  )  )
  
  print( summary( sim_list ) )
  
  return( sim_list )  
}

#' Show all mutations present in the database for a selected cancer cell line and reference Genome
#' @export
show_contained_mutations_for_cl = function( cl_name, ref_gen = "HG19"){
  
  suppressPackageStartupMessages(library("plyr"))
  suppressPackageStartupMessages(library("dplyr"))
  suppressPackageStartupMessages(library("stringr"))
  
  db_folder      = system.file("", package="Uniquorn")
  database_path  =  paste( db_folder, "uniquorn_db.sqlite3", sep ="/" )
  sim_list = as.data.frame( tbl( src_sqlite( database_path ), "sim_list_df" ), n = -1 )
  
  sim_list = sim_list[ sim_list$Ref_Gen == ref_gen,  ]
  mapping  = which( sim_list$CL %in% cl_name, arr.ind = T  )
  sim_list = sim_list[ mapping,  ]
  
  if ( length( mapping ) == 0  ){
    
    message(paste0("Could not find the cancer cell line ",cl_name, " in the database."), collapse= "")
    
  } else {
    
    print( paste0( c("Found ", dim(sim_list)[1], " many mutations for cancer cell line", cl_name  ," for reference genome ", ref_gen ), collapse = ""  )  )
    
    print( summary( sim_list ) )  
  }
  
  sim_list = sim_list[ mapping, ]
  return(sim_list)
}

#' Show all cancer cell lines in the database which contained the specified mutation and reference Genome. Closed interval coordinates. Format mutation: CHR_START_STOP, e.g. 1_123_123
#' @export
show_which_cls_contain_mutation = function( mutation_name, ref_gen = "HG19"){
  
  suppressPackageStartupMessages(library("plyr"))
  suppressPackageStartupMessages(library("dplyr"))
  suppressPackageStartupMessages(library("stringr"))
  
  db_folder      = system.file("", package="Uniquorn")
  database_path  =  paste( db_folder, "uniquorn_db.sqlite3", sep ="/" )
  sim_list = as.data.frame( tbl( src_sqlite( database_path ), "sim_list_df" ), n = -1 )
  
  sim_list = sim_list[ sim_list$Ref_Gen == ref_gen,  ]
  mapping  = which( sim_list$Fingerprint %in% mutation_name, arr.ind = T)
  sim_list = sim_list[ mapping,  ]
  
  if ( length( mapping ) == 0  ){
    
    message(paste0("Could not find any cancer cell line for the mutation ",mutation_name, " in the database."), collapse= "")
    
  } else {
    
    print( paste0( c("Found ", dim( sim_list )[1], " many cancer cell lines for mutation ", mutation_name  ," for reference genome ", ref_gen ), collapse = ""  )  )
    
    print(sim_list)
    print( summary( sim_list ) )  
  }
  
  return( sim_list )  
}