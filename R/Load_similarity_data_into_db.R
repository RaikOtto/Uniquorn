#' Loads similarity data into db
#' @export
load_similarity_data_into_db = function( similarity_matrix_data, db_path ){
  
  db_path = paste( db_path, "inst", sep = "/" )

  if ( ! dir.exists( db_path )  )
    dir.create( db_path )
  
  db_path  = paste( db_path,"Younikorn.db", sep ="/")
  print( paste0( "Storing data in db: ",db_path) )
  
  require( RSQLite )
  drv = dbDriver("SQLite")
  full_con = dbConnect( drv, dbname = db_path )
  
  dbWriteTable( 
    full_con,
    "similarity_matrix",
    as.data.frame( similarity_matrix_data ),
    overwrite = T
  )
}