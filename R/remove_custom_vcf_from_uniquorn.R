#' Removes a cancer cell line training fingerprint (vcf file) from the database. The names of all training sets can 
#' be seen by using the function \code{show_contained_cls}.
#' @param name_cl name of the cancer cell line training fingerprintt
#' @param ref_gen Reference genome version. All training sets are associated with a reference genome version. Default: GRCH37
#' @param distinct_mode Delete the CL sample from the database that is normalized together or separately. Choices: TRUE and FALSE
#' @import DBI
#' @usage 
#' remove_custom_vcf_from_database( 
#' 
#' name_cl, 
#' 
#' ref_gen = "GRCH37", 
#' 
#' distinct_mode = TRUE )
#' @return Message that indicates if the removal was succesful
#' @export
remove_custom_vcf_from_database = function( 
    name_cl,
    ref_gen = "GRCH37",
    distinct_mode = TRUE
){
    
    # pre processing

    name_cl = str_to_upper(name_cl)
    
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
    
    sim_list       = as.data.frame( DBI::dbReadTable( con, "sim_list") )
    sim_list_stats = as.data.frame( DBI::dbReadTable( con, "sim_list_stats") )
    dbDisconnect(con)
    
    if ( sum( grepl( name_cl, sim_list_stats$CL ) ) == 0 ){
        
        stop(paste0("No training set for a cancer cell line found for the name: ", name_cl))
        
    } else {
        
        print( 
            paste0(
                c("Found CL ",name_cl,
                ". Removing from database and recalculating training-sets."
            ), 
            collapse = "") 
        )
               
    }
        
    
    sim_list = sim_list[ sim_list$CL != name_cl,] # exclude sample here
    sim_list = sim_list[, which( colnames(sim_list) != "Ref_Gen"  ) ]
    sim_list = sim_list[, which( colnames(sim_list) != "Weight"  ) ]
    
    print("Found & removed the sample. Re-calculating the Cancer cell line data")
    
    list_of_cls = unique( as.character( sim_list$CL ) )
    panels = sapply( list_of_cls, FUN = stringr::str_split, "_"  )
    panels = as.character(unique( as.character( sapply( panels, FUN = utils::tail, 1) ) ))
    
    if (!distinct_mode){
        panels = paste0( c(panels), collapse ="|"  )
        database_path =  paste( package_path, "uniquorn_non_distinct_panels_db.sqlite3", sep ="/" )
    }
    
    print( paste( "Distinguishing between panels:",paste0( c(panels), collapse = ", "), sep = " ") )
    
    for (panel in panels) {
        
        print(panel)
        
        sim_list_panel   = sim_list[ grepl( panel, sim_list$CL) , ]
        member_var_panel = rep( 1, dim(sim_list_panel)[1] )
        
        sim_list_stats_panel = stats::aggregate( member_var_panel , by = list( sim_list_panel$CL ), FUN = sum )
        colnames(sim_list_stats_panel) = c( "CL", "Count" )
        
        print("Aggregating over mutational frequency to obtain mutational weight")
        
        weights_panel = stats::aggregate( member_var_panel , by = list( sim_list_panel$Fingerprint ), FUN = sum )
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
    
    drv = RSQLite::SQLite()
    con = DBI::dbConnect(drv, dbname = database_path)
    
    DBI::dbWriteTable( con, "sim_list", sim_list_global, overwrite = T )
    DBI::dbWriteTable( con, "sim_list_stats", sim_list_stats_global, overwrite = T )
    
    dbDisconnect(con)
    
    print(paste0( "Finished removing CL ", name_cl ))
}