#' Tutorial plots
#' @export
tutorial_analysis = function( db_path = system.file("database", package="Younikorn")  ){

  # init connection
  
  library("CancerCellLines")
  db_path  = paste0(db_path,"/Younikorn.db")
  print( paste0( "Parsing data and storing in db: ",db_path) )
  full_con = setupSQLite(db_path)
  dplyr_con = src_sqlite(full_con@dbname)
  
  # Example 1: Melanoma heatmap with MEK and BRAF inhibitors
  
  ex1_genes = c( 'BRAF', 'NRAS', 'CRAF', 'TP53')
  ex1_cell_lines <- dplyr_con %>% tbl('ccle_sampleinfo') %>% dplyr::filter(Site_primary=='skin') %>% collect %>% as.data.frame
  ex1_cell_lines <- ex1_cell_lines$CCLE_name
  print( ex1_cell_lines[1:10] )
  
  cat ("Press [enter] to continue")
  line <- readline()
  
  ex1_cls <- c( 'A2780')
  ex1_tall_df <- makeTallDataFrame(full_con, ex1_genes, ex1_cell_lines, ex1_drugs)
  
  print (ex1_tall_df)
}

