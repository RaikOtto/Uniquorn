#' Adds a custom vcf file to the three existing cancer cell line panels
#' It is strongly recommended to filter the SNPs with a minor allel frequency of more than 0.01.
#' @export
add_custom_vcf_to_database = function( 
    vcf_file_path,
    ref_gen = "GRCH37",
    name_cl = "",
    safe_mode = FALSE
    ){
    
    # pre processing
    require( "stringr", quietly = TRUE, warn.conflicts = FALSE )
    require( "RSQLite", quietly = TRUE, warn.conflicts = FALSE )
    require( "DBI",     quietly = TRUE, warn.conflicts = FALSE )
  
    package_path    = system.file("", package="Uniquorn")
    
    if (distinct_mode)
        database_path     =  paste( package_path, "uniquorn_distinct_panels_db.sqlite", sep ="/" )
    if (!distinct_mode)
        database_path     =  paste( package_path, "uniquorn_non_distinct_panels_db.sqlite", sep ="/" )
    
    # reading file
    vcf_fingerprint = parse_vcf_file( vcf_file )
    
    if ( output_file == ""  ){
        
        output_file = paste( vcf_file, "uniquorn_ident.tab", sep ="_")
        
    }else if ( dir.exists( output_file ) ){
        
        vcf_file_name = tail( as.character( unlist( str_split( vcf_file, "/" ) ) ),1 )
        output_file = paste(output_file,paste( vcf_file_name, "uniquorn_ident.tab", sep ="_"), sep = "/")
    }
    
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
    
    print("Finished reading database, identifying CL")
  
    if ( file.exists(vcf_file_path)  ){

        print( paste0( "Creating fingerprint from custom VCF file ", vcf_file_path  ) )
        
        if ( name_cl == "" ){
          
          vcf_identifier = as.character( tail( unlist( str_split( vcf_file_path, "/" ) ), 1) )
          name_cl = head( unlist( str_split( vcf_identifier, ".vcf|.VCF" ) )  , 1)
          name_cl = paste(name_cl, "CUSTOM", sep ="_"  )
          print( paste0( "No cl name provided, adding auto-generated fingerprint: ", name_cl ) )
          
        } else {
          
          name_cl = paste(name_cl, "custom", sep ="_"  )
          print( paste0( "Adding fingerprint with user-defined name: ", name_cl ) )
        }
        
        print( paste0( "Building fingerprint from file ",  vcf_file_path )  )
        vcf_fingerprint = as.character( parse_vcf_file( vcf_file_path ) )
        
        print(paste0("Adding fingerprint to reference genome: ", ref_gen))
        
        if (  sum( sim_list$CL == name_cl  ) != 0 ){
          
          stop( paste0 ( c( "Cell line ", name_cl,  " already present in dataset. Please remove CL from the dataset, using the respective function 'remove_cl_from_uniquorn_db' form the package or alter the name of the new CL."  ), collapse= ""  ) )
        } else {
            
          if ( safe_mode ){
            
            print( "Using safe mode" )
            
            vcf_fingerprint = vcf_fingerprint[ as.integer( which( vcf_fingerprint %in% sim_list$Fingerprint) ) ]
          }
          
          cl_vec          = rep( name_cl, length(vcf_fingerprint)  )
          weight_vec      = rep( 0.0,     length(vcf_fingerprint)  )
          ref_gen_vec     = rep( ref_gen, length(vcf_fingerprint)  )
          vcf_fingerprint = cbind( vcf_fingerprint, cl_vec, weight_vec, ref_gen_vec  )
          colnames(vcf_fingerprint ) = c("Fingerprint","CL","Weight","Ref_Gen")
          
          sim_list = rbind( sim_list, vcf_fingerprint )
          
          # calculate weights
          
          print( "Calculating weights for found mutations" )
          
          member_var = rep(1, dim(sim_list)[1] )
          weights = aggregate( member_var , by = list( sim_list$Fingerprint ), FUN = sum )
          weights$x = 1.0 / as.double( weights$x )
          
          mapping = match( as.character( sim_list$Fingerprint ), as.character( weights$Group.1) )
          sim_list$Weight = weights$x[mapping]
          
          print( "Calculating statistics for single CLs"  )
          
          sim_list_stats = aggregate( member_var , by = list( sim_list$CL ), FUN = sum )
          colnames(sim_list_stats) = c( "CL", "Count" )
          
          aggregation_all = stats::aggregate(
              x  = as.double( sim_list$Weight ),
              by = list( as.character( sim_list$CL ) ),
              FUN = sum
          )
          
          mapping_agg_stats = which( aggregation_all$Group.1 %in% sim_list_stats[,1], arr.ind = T  )
          sim_list_stats = cbind( sim_list_stats, aggregation_all$x[mapping_agg_stats] )
          
          Ref_Gen = rep( ref_gen, dim(sim_list_stats)[1]  )
          sim_list_stats = cbind( sim_list_stats, Ref_Gen )
          colnames( sim_list_stats ) = c( "CL","Count","All_weights","Ref_Gen" )
          
          ###
          
          print("Writing to DB")
          
          drv = RSQLite::SQLite()
          con = DBI::dbConnect( drv, dbname = database_path )
          
          DBI::dbWriteTable( con, "sim_list", sim_list, overwrite = T )
          DBI::dbWriteTable( con, "sim_list_stats", sim_list_stats, overwrite = T )
          dbDisconnect(con)
          
          ###
          
          print("Finished writing to database")
        }
    } else {
    
    message( paste0( "Did not find file: ", vcf_file_path)  )
    } 
    print("Finished")
}