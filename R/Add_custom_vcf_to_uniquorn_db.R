#' Adds a custom vcf file to the three existing cancer cell line panels
#' It is strongly recommended to filter the SNPs with a minor allel frequency of more than 0.01.
#' @export
add_custom_vcf_to_uniquorn_db = function( vcf_file_path, ref_gen = "hg19", name_cl = "" ){
    
  library("stringr")
  panels = c("CELLMINER","CCLE","COSMIC")
  ref_gen_path = paste( system.file("", package = "Uniquorn"), ref_gen, sep ="/" )
  
  if ( file.exists(vcf_file_path)  ){
  
    print( paste0( "Creating fingerprint from custom VCF file ", vcf_file_path  ) )
    
    if ( name_cl == "" ){
      
      vcf_identifier = as.character( tail( unlist( str_split( vcf_file_path, "/" ) ), 1) )
      name_cl = head( unlist( str_split( vcf_identifier, ".vcf|.VCF" ) )  , 1)
      name_cl = paste(name_cl, "custom", sep ="_"  )
      print( paste0( "No cl name provided, adding auto-generated fingerprint: ", name_cl ) )
      
    } else {
      
      print( paste0( "Adding fingerprint with user-defined name: ", name_cl ) )
    }
    
    print( paste0( "Building fingerprint from file ",  vcf_file_path )  )
    vcf_fingerprint = parse_vcf_file( vcf_file_path )
    
    print(paste0("Adding fingerprint to reference genome: ", ref_gen))
    
    for( panel in panels ){
      
      print( panel  )
      
      sim_list_file       = paste0( c( ref_gen_path, "/", "Fingerprint_",       panel, ".tab" ), collapse = "" )
      sim_list            = read.table( sim_list_file, sep = "\t", header = T)
      sim_list_stats      = read.table( sim_list_stats_file, sep = "\t", header = T)
      
      if ( ! name_cl %in% sim_list$CL   ){
      
        sim_list_stats_file = paste0( c( ref_gen_path, "/", "Fingerprint_stats_", panel, ".tab" ), collapse = "" )
        
        mapping_found_mutations = match( vcf_fingerprint, sim_list$Fingerprint, nomatch = 0  ) # update weights
        new_mutations       = mapping_found_mutations == 0
        
        add_data        <<- data.frame( 
          "Fingerprint" = vcf_fingerprint[new_mutations],
          "CL"          = rep( name_cl, length( vcf_fingerprint[new_mutations]  )  ),
          "Weight"      = rep(1 , length(vcf_fingerprint[new_mutations])  )
        )
        
        sub_matrix          = sim_list[ mapping_found_mutations, ]
        update_mutations    = unique( as.character( sub_matrix$Fingerprint ) )
        
        update_weights = function(  single_fingerprint, sub_matrix, add_data ){
          
          mapping = which( as.character( sub_matrix$Fingerprint ) == single_fingerprint  )
          weight  = 1.0 / ( length( mapping ) + 1 )
          
          sub_matrix$Weight[mapping] = weight
          
          add_data = data.frame(
            "Fingerprint" = c( as.character( add_data$Fingerprint ), as.character( single_fingerprint ) ),
            "CL"          = c( as.character( add_data$CL ), name_cl ),
            "Weight"      = c( add_data$Weight, as.character( weight ) )
          )
          
        }
        
        lapply( update_mutations, FUN = update_weights, sub_matrix, add_data )
        
        stats_update = matrix( ncol = 2, c(name_cl, length( vcf_fingerprint )  )  )
        colnames(stats_update) = c( "CL", "Count" )
        sim_list_stats = rbind( sim_list_stats, stats_update )
        
        sim_list$Weight[ which( sim_list$Fingerprint %in% sub_matrix$Fingerprint )  ] = as.character( sub_matrix$Weight )
        
        sim_list = rbind( sim_list, add_data )
        
        write.table( sim_list, file = sim_list_file, sep ="\t", quote= F, row.names = F  )
        write.table( sim_list_stats, file =  sim_list_stats_file, sep ="\t", quote= F, row.names =F  )
      
      } else {
        
        message( paste0 ( c( "Cell line ", name_cl,  " already present in dataset. Please remove CL from the dataset, using the respective function 'remove_cl_from_uniquorn_db' form the package or alter the name of the new CL."  ), collapse= ""  ) )
      }
    }
  } else {
    
    message( paste0( "Did not find file: ", vcf_file_path)  )
  } 
  print("Finished")
}