#' Create similarity matrix from raw_data
#' @export
create_similarity_matrix = function( raw_data  ){
  
  library("stringr")
  
  unique_list_cls = unique( raw_data$CL_ident  )
  n_cols = length( unique_list_cls )
  
  coords = apply( 
    as.matrix( 
      data.frame( 
        raw_data$Chr,
        raw_data$start,
        raw_data$stop
      )
    ),
    MARGIN = 1,
    FUN = paste0,
    collapse = "_"
  )
  
  unique_list_coords = unique( coords )
  n_rows = length( unique_list_coords  )
  
  res = matrix( rep( "0", n_rows * n_cols) , nrow = n_rows, ncol = n_cols  )
  
  return_membership_status_mutation_cl = function( mutation ){
    
    cls_matching = unique( raw_data$CL_ident[ coords == mutation ]  )
    res[ which( unique_list_coords == mutation ) ,  which( unique_list_cls %in% cls_matching )  ] <<- "1"
  }
  
  lapply( unique_list_coords, FUN = return_membership_status_mutation_cl )

  colnames(res) = unique_list_cls
  rownames(res) = unique_list_coords
  
  return( res )
}