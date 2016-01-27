#' inititate_db_and_load_data
#' 
#' Intern utility function, loads database and return the sim_list and sim_list_stats variables.
#' 
#' @inheritParams identify_vcf_file
#' @param request_tables Names of the tables to be extracted from the database
#' @return the sim_list and sim_list_stats variable
#' @usage 
#' inititate_db_and_load_data( ref_gen = "GRCH37", distinct_mode = TRUE )
#' @import DBI RSQLite
inititate_db_and_load_data = function( ref_gen, distinct_mode, request_tables ){
    
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
    
    res = c()
    
    for ( request_table in request_tables ){
        res = c(res, as.data.frame( DBI::dbReadTable( con, request_table) ) )
    }
    DBI::dbDisconnect(con)
    
    return( res )
}

#' write_data_to_db
#' 
#' Intern utility function, writes to database the sim_list and sim_list_stats variables
#' 
#' @param sim_list R Table which contain a mapping of mutations to cancer cell lines for a specific reference genome
#' @param sim_list_stats Contains an aggergated R table that shows the amount and weight of mutations for a reference genome
#' @return the sim_list and sim_list_stats variable
#' @usage 
#' write_data_to_db( sim_list = sim_list, sim_list_stats sim_list_stats, ref_gen = "GRCH37", distinct_mode = TRUE )
#' @import DBI RSQLite
write_data_to_db = function( sim_list, sim_list_stats, ref_gen, distinct_mode ){
    
    package_path    = system.file("", package="Uniquorn")
    
    if (distinct_mode)
        database_path =  paste( package_path, "uniquorn_distinct_panels_db.sqlite", sep ="/" )
    if (!distinct_mode)
        database_path =  paste( package_path, "uniquorn_non_distinct_panels_db.sqlite", sep ="/" )
    
    if( ! file.exists( database_path ) ){
        warning("Writing to default database! This is not recommended. Please consider adding the CCLE and CoSMIC training sets for optimal Uniquorn function.")
        database_default_path =  paste( package_path, "uniquorn_db_default.sqlite", sep ="/" )
        database_path = database_default_path
    }
    
    drv = RSQLite::SQLite()
    con = DBI::dbConnect(drv, dbname = database_path)
    
    
    DBI::dbDisconnect(con)

}