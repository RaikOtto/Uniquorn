#' initiate_db_and_load_data
#' 
#' Intern utility function, loads database and return the sim_list and sim_list_stats variables.
#' 
#' @inheritParams identify_vcf_file
#' @param request_table Names of the tables to be extracted from the database
#' @param load_default_db Indicate whether the default db should be used as source for the data
#' @return Returns the sim_list and sim_list_stats variable
#' @usage 
#' initiate_db_and_load_data(
#' ref_gen, 
#' request_table,
#' load_default_db )
#' @import DBI RSQLite
initiate_db_and_load_data = function( ref_gen, request_table, load_default_db = FALSE ){
    
    package_path = system.file("", package="Uniquorn")
    default_database_path =  paste( package_path, "uniquorn_db_default.sqlite", sep ="/" )

    database_path = paste( package_path, "uniquorn_distinct_panels_db.sqlite", sep ="/" )

    drv = RSQLite::SQLite()
    
    if( file.exists( database_path) & (! load_default_db ) & (file.size( database_path ) != 0 ) )
        con = DBI::dbConnect(drv, dbname = database_path)
    else
        con = DBI::dbConnect(drv, dbname = default_database_path)

    res = as.data.frame( DBI::dbReadTable( con, request_table) )
    
    DBI::dbDisconnect(con)
    
    return( res )
}

#' write_data_to_db
#' 
#' Intern utility function, writes to database the sim_list and sim_list_stats variables
#' 
#' @inheritParams identify_vcf_file
#' @param content_table Tables to be written in db
#' @param table_name Name of the table to be written into the DB
#' @param overwrite Overwrite the potentially existing table
#' @param test_mode Is this a test? Just for internal use 
#' @return the sim_list and sim_list_stats variable
#' @usage 
#' write_data_to_db( 
#' content_table, 
#' table_name,
#' ref_gen,
#' overwrite,
#' test_mode )
#' @import DBI RSQLite
write_data_to_db = function( 
    content_table, 
    table_name, 
    ref_gen = "GRCH37",
    overwrite = TRUE, 
    test_mode = FALSE 
    ){
    
    package_path    = system.file("", package="Uniquorn")
    
    database_path =  paste( package_path, "uniquorn_distinct_panels_db.sqlite", sep ="/" )

    drv = RSQLite::SQLite()
    con = DBI::dbConnect(drv, dbname = database_path)
    
    if( test_mode )
        con = DBI::dbConnect(drv, dbname = "")
    
    DBI::dbWriteTable( con, name = table_name, value = as.data.frame( content_table ), overwrite = overwrite )
    
    DBI::dbDisconnect(con)

}

#' Re-calculate sim_list_weights
#' 
#' This function re-calculates the weights of mutation after a change of the training set
#' @inheritParams write_data_to_db
#' @return A list containing both the sim_list at pos 1 and sim_list_stats at pos 2 data frames.
#' @param sim_list R Table which contains a mapping 
#' from mutations/ variations to their containing CLs
re_calculate_cl_weights = function( sim_list, ref_gen ){
    
    package_path    = system.file("", package="Uniquorn")
    
    list_of_cls = unique( as.character( sim_list$CL ) )
    panels = sapply( list_of_cls, FUN = stringr::str_split, "_"  )
    panels = as.character(unique( as.character( sapply( panels, FUN = utils::tail, 1) ) ))
    
    print( paste( "Distinguishing between panels:",
        paste0( c(panels), collapse = ", "), sep = " ") )
    
    if(! exists("sim_list_global"))
        sim_list_global = sim_list[0,]
    
    for (panel in panels) {
        
        print(panel)
        
        sim_list_panel   = sim_list[ grepl( panel, sim_list$CL) , ]
        member_var_panel = rep( 1, dim(sim_list_panel)[1] )
        
        sim_list_stats_panel = stats::aggregate( member_var_panel , 
            by = list( sim_list_panel$CL ), FUN = sum )
        
        colnames(sim_list_stats_panel) = c( "CL", "Count" )
        
        print("Aggregating over mutational frequency to obtain mutational weight")
        
        weights_panel = stats::aggregate( member_var_panel , 
            by = list( sim_list_panel$Fingerprint ), FUN = sum )
        
        weights_panel$x = 1.0 / as.double( weights_panel$x )
        
        mapping_panel = match( as.character( sim_list_panel$Fingerprint ), 
            as.character( weights_panel$Group.1) )
        
        sim_list_panel = cbind( sim_list_panel, weights_panel$x[mapping_panel] )
        colnames( sim_list_panel )[3] = "Weight"
        
        # calculate weights
        
        aggregation_all_panel = stats::aggregate( 
            x  = as.double( sim_list_panel$Weight ),
            by = list( as.character( sim_list_panel$CL ) ),
            FUN = sum
        )
        
        mapping_agg_stats_panel = which( aggregation_all_panel$Group.1 %in% 
            sim_list_stats_panel[,1], arr.ind = TRUE  )
        sim_list_stats_panel = cbind( sim_list_stats_panel, 
            aggregation_all_panel$x[mapping_agg_stats_panel] )
        
        #print("Finished aggregating, writing to database")
        
        Ref_Gen = rep(ref_gen, dim(sim_list_panel)[1]  )
        sim_list_panel = cbind( sim_list_panel, Ref_Gen )
        Ref_Gen = rep( ref_gen, dim(sim_list_stats_panel)[1]  )
        sim_list_stats_panel = cbind( sim_list_stats_panel, Ref_Gen )
        colnames( sim_list_stats_panel ) = c( "CL","Count","All_weights","Ref_Gen" )
        
        sim_list_global = rbind(sim_list_global,sim_list_panel)
        
        if(! exists("sim_list_stats_global"))
            sim_list_stats_global = sim_list_stats_panel[0,]
        
        sim_list_stats_global = rbind( sim_list_stats_global, sim_list_stats_panel  )
        
    }
    
    return( list( as.data.frame( sim_list_global ), as.data.frame( sim_list_stats_global ) ) )
}