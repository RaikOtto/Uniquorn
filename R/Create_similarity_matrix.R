#' Create similarity matrix from raw_data
#' @export
create_similarity_matrix = function( fingerprint_data, cl_data ){
  
  library("stringr")
  
  nr_fingerprints = dim(fingerprint_data)[1]
  nr_cls          = dim(cl_data)[1]
  
  adv = 0

  create_fingerprint = function( cl_entry ,fp_list  ){
   
    fp = unlist( 
      str_split( 
        as.character( 
          cl_entry
        ),
        pattern = ","
      )
    )
    
    res = match( 
      as.vector(
        fp_list),
      fp,
      nomatch = 0
    )
    res[ res != 0  ] = 1 
    
    
    stat = round( (adv / as.double(nr_cls)) * 100, 1 )
    adv <<- adv + 1
      
    if ( stat != round( (adv / as.double(nr_cls)) * 100, 1 ) )
      print( paste( round( (adv / as.double(nr_cls)) * 100, 1 ), "% finished", sep =" " ) )

    return(res)
  }
  
  sim_list = lapply(
    cl_data$Fingerprints,
    FUN = create_fingerprint,
    fingerprint_data$Fingerprint
  )
  
  #res_common = matrix( as.integer( unlist( sim_list ) ), ncol = length( sim_list)  )
  rownames( res_common ) = fingerprint_data$Fingerprint
  colnames( res_common ) = cl_data$CL
  
  return( sim_list )
}