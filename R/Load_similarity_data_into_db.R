#' Loads similarity data into db
load_similarity_data_into_db = function( similarity_matrix_data, db_path = system.file("", package="Younikorn") ){
  
  db_path = paste( db_path, "inst", sep = "/" )

  if ( ! dir.exists( db_path )  )
    dir.create( db_path )
  
  db_path  = paste( db_path,"Younikorn.db", sep ="/")
  print( paste0( "Storing data in db: ",db_path) )
  
  require( RSQLite )
  suppressMessages( require(dplyr  ) )
  
  file.remove(db_path)
  similarity_db = src_sqlite( db_path, create = T)
  
  similarity_matrix = as.data.frame( similarity_matrix_data )
  
  full_con = copy_to(
    similarity_db, 
    df = similarity_matrix, 
    temporary = F,
    
    indexes = list( 
      "mutational_similarity_marker"
    )
  )
}