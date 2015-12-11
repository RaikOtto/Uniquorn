#' Parses the vcf file and predicts the identity of the sample
#' @export
identify_vcf_file = function( vcf_file_path, output_path = "", panels = c("CELLMINER","CCLE","COSMIC"), type = "unique" ){
  # types = c("unique","non_unique")
  library("stringr")
  
  print( paste0( "Creating fingerprint from VCF file ", vcf_file_path  ) )
  vcf_fingerprint = parse_vcf_file( vcf_file_path )
  
  sim_list_store_mat = matrix( list(), ncol = length( type ), nrow = length( panels ) , dimnames = list( panels, type ) )
  
  if ( output_path == ""  ){
    output_file = paste( vcf_file_path, "ident.tab", sep ="_")
  } else {
    output_file = paste( paste0( output_path, tail( unlist(str_split(vcf_file_path, "/" ) ), 1 )  ), ".tab", sep ="_ident")
  }
  
  res_table <<- data.frame(
    
    "CL" = as.character(),
    "Found_muts_abs" = as.character(),
    "Found_muts_rel" = as.character(),
    "Count_mutations_cl" = as.character(),
    "Passed_threshold" = as.character()
  )

  for( panel in panels ){
    
    sim_list_file = paste( system.file("", package = "Uniquorn"), 
      paste0(
      c(
        type,
       "_parsed_DB_mut_labels_",
       panel,
       ".tab"
      ),
        collapse= ""
      ),
      sep = "/"
    )

    print( paste0( "Loading similarity data from file ",  sim_list_file )  )
    
    if ( length( sim_list_store_mat[[ panel, type ]])  != 0 ){
        
      sim_list = sim_list_store_mat[[ panel, type ]]
      
    } else {
      
      if (type == "non_unique"){
        
        # pass
        
      } else {
        
        sim_list = read.table( sim_list_file, sep = "\t", header = T)

      }
      sim_list_store_mat[[ panel, type ]] = sim_list

    }

    print( panel  )
    
    mapping = which( sim_list$Fingerprint %in% vcf_fingerprint )
    
    candidate_hits_abs = table( sim_list$CL[mapping] )
    cl_match = match( names(table( sim_list$CL[mapping] )), sim_list$CL )
    candidate_hits_rel = round( candidate_hits_abs / sim_list$Count[cl_match], 3 ) * 100
    
    nr_cls = length( unique( sim_list$CL  )  )
    passed_threshold_vec = rep( F, nr_cls )
    passed_threshold_vec[ ( candidate_hits_abs >= 2 ) & ( candidate_hits_rel >= 3 ) ] = T
    
    #"Found_muts_weighted" = as.double( res_intersect[2,] ),
    #"Found_muts_rel_weighted" = as.double( res_intersect[2,] ),
    #"Passed_threshold_weighted" = rep( F, nr_cls )

    # don't sort me bro!
    
    res_table = data.frame( 
      "CL"                 = c( as.character( res_table$CL) , as.character( sim_list$CL[cl_match] ) ),
      "Found_muts_abs"     = c( as.character( res_table$Found_muts_abs), as.character( candidate_hits_abs ) ),
      "Found_muts_rel"     = c( as.character( res_table$Found_muts_rel), as.character(  candidate_hits_rel ) ),
      "Count_mutations_cl" = c( as.character( res_table$Count_mutations_cl), as.character(  sim_list$Count[cl_match] ) ),
      "Passed_threshold"   = c( as.character( res_table$Passed_threshold), as.character( passed_threshold_vec ) )
    )

  }
  
  res_table = res_table[ order( as.double( res_table$Found_muts_rel ), decreasing = T),  ]
  
  if ( dim( res_table[ res_table$Passed_threshold ,])[1] == 1 ){
    
    print( paste0( "Candidate(s): ", paste0( ( res_table$CL_name[ res_table$Passed_treshold  ] ), collapse = "," ) )  )
    
  } else if ( dim( res_table[ res_table$Passed_threshold ,])[1] > 1 ){
    
    print( paste0( "Candidates: ", paste0( ( res_table$CL_name[ res_table$Passed_treshold  ] )[1:2], collapse = "," ) )  )
    
  } else {
    
    print( paste0("No CL mutational fingerprint with sufficient similarity found." ) )
  }  
  
  print( paste0("Storing information in table: ",output_file ) )
  write.table( res_table, output_file, sep ="\t", row.names = F, quote = F  )
}