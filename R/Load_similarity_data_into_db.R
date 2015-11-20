#' Loads similarity data into db
#' @export
load_similarity_data_into_db = function( sim_list, fingerprint_data, db_path){

  print( paste0( "Storing data in db: ",db_path) )
  
  require( RSQLite )
  suppressMessages( require( dplyr ) )
  
  if ( file.exists( db_path ) )
    file.remove(db_path)
  
  similarity_db = src_sqlite( db_path, create = T)
  
  similarity_matrix = as.data.frame( 
    matrix( unlist( sim_list ) , ncol = length(sim_list) )
  )
  similarity_matrix = rbind( similarity_matrix, fingerprint_data$Fingerprint )

  colnames( similarity_matrix )[1] = "Mutational_fingerprint"
  rownames( similarity_matrix ) = fingerprint_data$Fingerprint
  
  copy_to(
    similarity_db, 
    df = similarity_matrix, 
    temporary = F,
    
    indexes = list(
      colnames( similarity_matrix_data )
    )
  )
}

load_data_into_db = function( fingerprint_data, cl_data, db_path){
  
  print( paste0( "Storing data in db: ",db_path) )
  
  require( RSQLite )
  suppressMessages( require( dplyr ) )
  
  if ( file.exists( db_path ) )
    file.remove(db_path)
  
  similarity_db = src_sqlite( db_path, create = T)
  
  similarity_matrix = as.data.frame( 
    matrix( unlist( sim_list ),
    ncol = length(sim_list) )
  )
  
  similarity_matrix = rbind( similarity_matrix, fingerprint_data$Fingerprint )
  
  colnames( similarity_matrix )[1] = "Mutational_fingerprint"
  rownames( similarity_matrix ) = fingerprint_data$Fingerprint
  
  copy_to(
    similarity_db, 
    df = similarity_matrix, 
    temporary = F,
    
    indexes = list(
      colnames( similarity_matrix_data )
    )
  )
}