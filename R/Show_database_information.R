#' Show all cancer cell line identifier present in the database for a selected reference genome.
#' @param ref_gen Reference genome version. All training sets are associated with a reference genome version. Default: GRCH37
#' @param distinct_mode Show training data for the commonly or separately normalized training sets. Options: TRUE/ FALSE
#' @return R table which contains the identifier of all cancer cell line samples with the specific reference genome
#' @import DBI stringr RSQLite
#' @export
show_contained_cls = function( ref_gen = "GRCH37", distinct_mode = T ){

    print(paste0("Reference genome: ",ref_gen))
    
    package_path    = system.file("", package="Uniquorn")
    
    if (distinct_mode)
        database_path =  paste( package_path, "uniquorn_distinct_panels_db.sqlite", sep ="/" )
    if (!distinct_mode)
        database_path =  paste( package_path, "uniquorn_non_distinct_panels_db.sqlite", sep ="/" )
    
    if( ! file.exists( database_path ) ){
        warning("Only found the vanilla CellMiner default database")
        database_default_path =  paste( package_path, "uniquorn_db_default.sqlite", sep ="/" )
        database_path = database_default_path
    }
        
    
    drv = RSQLite::SQLite()
    con = DBI::dbConnect(drv, dbname = database_path)

    sim_list_stats = as.data.frame( DBI::dbReadTable( con, "sim_list_stats") )
    dbDisconnect(con)
    
    print( paste0( c("Found ", dim(sim_list_stats)[1], " many cancer cell lines fingerprints for reference genome ", ref_gen ), collapse = ""  )  )

    print( paste( "CoSMIC CLP: ", as.character( sum( grepl( "_COSMIC", sim_list_stats$CL ) ) ) ) )
    print( paste( "CCLE: ", as.character( sum( grepl( "_CCLE", sim_list_stats$CL ) ) ) ) )
    print( paste( "CellMiner: ", as.character( sum( grepl( "_CELLMINER", sim_list_stats$CL ) ) ) ) )
    print( paste( "CUSTOM: ", as.character( sum( grepl( "_CUSTOM", sim_list_stats$CL ) ) ) ) )
    
    return( sim_list_stats )  
}

#' Show all mutations present in the database for a selected reference Genome
#' @param ref_gen Reference genome version
#' @param distinct_mode Show mutations for either distinct or non-distinct normalization of mutational weights
#' @return R Table which contains all mutations associated with a particular cancer cell line for a specified reference genome
#' @export
show_contained_mutations = function( ref_gen = "GRCH37", distinct_mode = T ){
  
    require("DBI", quietly = TRUE, warn.conflicts = FALSE)
    require("stringr", quietly = TRUE, warn.conflicts = FALSE)
    
    print(paste0("Reference genome: ",ref_gen))
    
    package_path    = system.file("", package="Uniquorn")
    
    if (distinct_mode)
        database_path =  paste( package_path, "uniquorn_distinct_panels_db.sqlite", sep ="/" )
    if (!distinct_mode)
        database_path =  paste( package_path, "uniquorn_non_distinct_panels_db.sqlite", sep ="/" )
    
    if( ! file.exists( database_path ) )
        stop(paste0("Did not find the database from which to delete the dataset: ", database_path))
    
    drv = RSQLite::SQLite()
    con = DBI::dbConnect(drv, dbname = database_path)
    
    sim_list = as.data.frame( DBI::dbReadTable( con, "sim_list") )
    dbDisconnect(con)
    
    print( paste0( c("Found ", dim(sim_list)[1], " many cancer cell lines associated mutations for reference genome ", ref_gen ), collapse = ""  )  )
  
    print( summary( sim_list ) )
  
    return( sim_list )  
}

#' Show all mutations present in the database for a selected cancer cell line and reference Genome
#' @param ref_gen Reference genome version
#' @param distinct_mode Show mutations for either distinct or non-distinct normalization of mutational weights
#' @param name_cl Name of the cancer cell line sample stored in the database
#' @import DBI stringr
#' @return R table which contains all mutations associated with the defined cancer cell line and reference genome
#' @export
show_contained_mutations_for_cl = function( name_cl, ref_gen = "GRCH37", distinct_mode = T){

    print(paste0("Reference genome: ",ref_gen))
    
    package_path    = system.file("", package="Uniquorn")
    
    if (distinct_mode)
        database_path =  paste( package_path, "uniquorn_distinct_panels_db.sqlite", sep ="/" )
    if (!distinct_mode)
        database_path =  paste( package_path, "uniquorn_non_distinct_panels_db.sqlite", sep ="/" )
    
    if( ! file.exists( database_path ) )
        stop(paste0("Did not find the database from which to delete the dataset: ", database_path))
    
    drv = RSQLite::SQLite()
    con = DBI::dbConnect(drv, dbname = database_path)
    
    sim_list = as.data.frame( DBI::dbReadTable( con, "sim_list") )
    dbDisconnect(con)
  
    sim_list = sim_list[ sim_list$Ref_Gen == ref_gen,  ]
    mapping  = which( sim_list$CL %in% name_cl, arr.ind = T  )
    sim_list = sim_list[ mapping,  ]
    
    if ( length( mapping ) == 0  ){
    
        message(paste0("Could not find the cancer cell line ",name_cl, " in the database."), collapse= "")
    
    } else {
    
        print( paste0( c("Found ", dim(sim_list)[1], " many mutations for cancer cell line", name_cl  ," for reference genome ", ref_gen ), collapse = ""  )  )
    }
    
    return(sim_list)
}

#' Show all cancer cell lines in the database which contained the specified mutation and reference Genome. Closed interval coordinates. Format mutation: CHR_START_STOP, e.g. 1_123_123
#' @param mutation_name Name of the mutation in the format CHROMOSOME_START_STOP, e.g. '11_244501_244510'
#' @param ref_gen Reference genome version
#' @param distinct_mode Show mutations for either distinct or non-distinct normalization of mutational weights
#' @import DBI stringr
#' @return R table which contains all cancer cell line samples which contain the specified mutation with respect to the specified reference genome version
#' @export
show_which_cls_contain_mutation = function( mutation_name, ref_gen = "GRCH37", distinct_mode = TRUE){
  
    print(paste0("Reference genome: ",ref_gen))
    
    package_path    = system.file("", package="Uniquorn")
    
    if (distinct_mode)
        database_path =  paste( package_path, "uniquorn_distinct_panels_db.sqlite", sep ="/" )
    if (!distinct_mode)
        database_path =  paste( package_path, "uniquorn_non_distinct_panels_db.sqlite", sep ="/" )
    
    if( ! file.exists( database_path ) )
        stop(paste0("Did not find the database from which to delete the dataset: ", database_path))
    
    drv = RSQLite::SQLite()
    con = DBI::dbConnect(drv, dbname = database_path)
    
    sim_list = as.data.frame( DBI::dbReadTable( con, "sim_list") )
    dbDisconnect(con)
  
    sim_list = sim_list[ sim_list$Ref_Gen == ref_gen,  ]
    mapping  = which( sim_list$Fingerprint %in% mutation_name, arr.ind = T)
    sim_list = sim_list[ mapping,  ]
    
    if ( length( mapping ) == 0  ){
    
        message(paste0("Could not find any cancer cell line for the mutation ",mutation_name, " in the database."), collapse= "")
    
    } else {
    
        print( paste0( c("Found ", dim( sim_list )[1], " many cancer cell lines for mutation ", mutation_name  ," for reference genome ", ref_gen ), collapse = ""  )  )
    }
    
    return( sim_list )  
}