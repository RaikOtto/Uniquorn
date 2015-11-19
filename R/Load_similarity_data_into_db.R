#' Loads similarity data into db
load_similarity_data_into_db = function( similarity_matrix_data, db_path = system.file("", package="Younikorn") ){

  print( paste0( "Storing data in db: ",db_path) )
  
  require( RSQLite )
  suppressMessages( require(dplyr  ) )
  
  file.remove(db_path)
  similarity_db = src_sqlite( db_path, create = T)
  
  similarity_matrix = as.data.frame( similarity_matrix_data )
  
  copy_to(
    similarity_db, 
    df = similarity_matrix, 
    temporary = F,
    
    indexes = list( 
      colnames( similarity_matrix_data )
    )
  )
}