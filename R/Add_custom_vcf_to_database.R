#' Adds a custom vcf file to the three existing cancer cell line panels
#' It is strongly recommended to filter the SNPs with a minor allel frequency of more than 0.01.
#' @export
add_custom_vcf_to_database = function( vcf_file_path, ref_gen = "GRCH37", name_cl = "", safe_mode = F ){
    
  # pre processing
  suppressPackageStartupMessages(library("plyr"))
  suppressPackageStartupMessages(library("dplyr"))
  suppressPackageStartupMessages(library("stringr"))
  
  db_folder    = system.file("", package="Uniquorn")
  database_path=  paste( db_folder, "uniquorn_db.sqlite3", sep ="/" )
  sim_list     = as.data.frame( tbl( src_sqlite( database_path ), "sim_list_df" ), n = -1 )
  
  if ( file.exists(vcf_file_path)  ){
  
    print( paste0( "Creating fingerprint from custom VCF file ", vcf_file_path  ) )
    
    if ( name_cl == "" ){
      
      vcf_identifier = as.character( tail( unlist( str_split( vcf_file_path, "/" ) ), 1) )
      name_cl = head( unlist( str_split( vcf_identifier, ".vcf|.VCF" ) )  , 1)
      name_cl = paste(name_cl, "custom", sep ="_"  )
      print( paste0( "No cl name provided, adding auto-generated fingerprint: ", name_cl ) )
      
    } else {
      
      name_cl = paste(name_cl, "custom", sep ="_"  )
      print( paste0( "Adding fingerprint with user-defined name: ", name_cl ) )
    }
    
    print( paste0( "Building fingerprint from file ",  vcf_file_path )  )
    vcf_fingerprint = as.character( parse_vcf_file( vcf_file_path ) )
    
    print(paste0("Adding fingerprint to reference genome: ", ref_gen))
    
    if (  sum( sim_list$CL == name_cl  ) != 0 ){
      
      message( paste0 ( c( "Cell line ", name_cl,  " already present in dataset. Please remove CL from the dataset, using the respective function 'remove_cl_from_uniquorn_db' form the package or alter the name of the new CL."  ), collapse= ""  ) )
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
      
      print("Writing to DB")
      
      if (file.exists(database_path))
        file.remove( database_path )
      uni_db   = src_sqlite( database_path, create = T )
      sim_list_df       = tbl_df( sim_list )
      sim_list_stats_df = tbl_df( sim_list_stats )
      
      copy_to( uni_db, sim_list_df, temporary = F, 
        indexes = list(
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
      print("Finished writing to database")
    }
  } else {
    
    message( paste0( "Did not find file: ", vcf_file_path)  )
  } 
  print("Finished")
}