#' Adds a custom vcf file to the three existing cancer cell line panels
#' @param vcf_file Input vcf file. Only one sample column allowed.
#' @param ref_gen Reference genome version. All training sets are associated with a reference genome version. Default: GRCH37
#' @param name_cl Name of the to-be-added cancer cell line sample. '_CUSTOM' will automatically be added as suffix.
#' @param safe_mode Only add mutations to the database where there already are mutations found in the cannonical cancer cell lines. This is a safety mechanism against overfitting if there are too few custom training samples.
#' @param distinct_mode Show training data for the commonly or separately normalized training sets. Options: TRUE/ FALSE
#' @import DBI stringr 
#' @export
add_custom_vcf_to_database = function( 
    vcf_file_path,
    ref_gen = "GRCH37",
    name_cl = "",
    safe_mode = FALSE,
    distinct_mode = TRUE
    ){
    
    # pre processing
    require( "stringr", quietly = TRUE, warn.conflicts = FALSE )
    require( "RSQLite", quietly = TRUE, warn.conflicts = FALSE )
    require( "DBI",     quietly = TRUE, warn.conflicts = FALSE )
    
    name_cl = str_to_upper(name_cl)
    
    print(paste0("Reference genome: ",ref_gen))
  
    package_path    = system.file("", package="Uniquorn")
    
    if (distinct_mode)
        database_path     =  paste( package_path, "uniquorn_distinct_panels_db.sqlite", sep ="/" )

    if (!distinct_mode)
        database_path     =  paste( package_path, "uniquorn_non_distinct_panels_db.sqlite", sep ="/" )
    
    if( ! file.exists( database_path ) ){
        
        database_path = paste( package_path, "uniquorn_db_default.sqlite", sep ="/" )
        warning("CCLE & CoSMIC CLP cancer cell line fingerprint NOT found, defaulting to 60 CellMiner cancer cell lines! 
                It is strongly advised to add ~1900 CCLE & CoSMIC CLs, see readme.")
    }
    
    drv = RSQLite::SQLite()
    con = DBI::dbConnect(drv, dbname = database_path)
    
    sim_list       = as.data.frame( DBI::dbReadTable( con, "sim_list") )
    sim_list_stats = as.data.frame( DBI::dbReadTable( con, "sim_list_stats") )
    
    sim_list = sim_list[ sim_list$Ref_Gen == ref_gen,]
    sim_list_stats = sim_list_stats[ sim_list_stats$Ref_Gen == ref_gen,]
    
    if ( dim(sim_list)[1] == 0 )
        warning( paste( "Warning: Identification might be spurious due to low amount of training sample. No cancer cell line data stored at this point for reference genome:", ref_gen, sep =" " ) )
    
    sim_list = sim_list[, which( colnames(sim_list) != "Ref_Gen"  ) ]
    sim_list = sim_list[, which( colnames(sim_list) != "Weight"  ) ]
    
    # reading vcf
    
    if ( name_cl == "" ){
        
        vcf_identifier = as.character( tail( unlist( str_split( vcf_file_path, "/" ) ), 1) )
        name_cl = head( unlist( str_split( vcf_identifier, ".vcf|.VCF" ) )  , 1)
        name_cl = str_to_upper( paste(name_cl, "CUSTOM", sep ="_"  ) )
        print( paste0( "No cl name provided, adding auto-generated fingerprint: ", name_cl ) )
        
    } else {
        
        name_cl = str_to_upper( paste(name_cl, "custom", sep ="_"  ) )
        print( paste0( "Adding fingerprint with user-defined name: ", name_cl ) )
    }
    
    print( paste0( "Building fingerprint from file ",  vcf_file_path )  )
    
    vcf_fingerprint = as.character( parse_vcf_file( vcf_file_path ) )
    
    if(safe_mode)
        vcf_fingerprint = vcf_fingerprint[ which( vcf_fingerprint %in% sim_list$Fingerprint ) ]
    
    vcf_fingerprint = data.frame( 
        "Fingerprint" = vcf_fingerprint,
        "CL"  = rep( name_cl, length(as.character(vcf_fingerprint)) )
    )
    
    sim_list = rbind( sim_list, vcf_fingerprint)
    
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
    
    drv = RSQLite::SQLite()
    con = DBI::dbConnect(drv, dbname = database_path)
    
    DBI::dbWriteTable( con, "sim_list", sim_list_global, overwrite = T )
    DBI::dbWriteTable( con, "sim_list_stats", sim_list_stats_global, overwrite = T )
    dbDisconnect(con)
    
    print("Finished adding CL")
}