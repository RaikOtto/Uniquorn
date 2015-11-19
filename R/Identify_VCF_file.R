#' Parses the vcf file and predicts the identity of the sample
#' @export
identify_vcf_file = function( vcf_file_path, db_path = system.file("", package="Younikorn" ) ){
  
  library( "dplyr" )
  library( "RSQLite" )
  
  source( "./R/Parse_VCF_data.R" )
  vcf_fingerprint = parse_vcf_file( vcf_file_path )
  
  #fres = data.frame( dplyr::filter( full_con, mutational_similarity_marker %in% c("1_10521399_10521399","1_10521520_10521520") ))
  #res2 = dplyr::arrange( full_con, "select * from similarity_matrix_data" )
  
  init_connection = function ( db_path ){

    if ( grepl( "Younikorn.db", c(db_path)) != T )
      db_path = paste( db_path, "inst/Younikorn.db", sep = "/" )

    drv = dbDriver("SQLite")
    con = dbConnect(
      drv,
      dbname = db_path
    )
    con = src_sqlite( con@dbname )
      
    return( con )
  }
  
  full_con = init_connection( db_path )
  
  raw_res = as.data.frame(

     dplyr::filter( 
      tbl( full_con, from = "similarity_matrix" ),
      mutational_similarity_marker %in% vcf_fingerprint
     )
  )
  cl_names = colnames(raw_res)[-1]

  res_common = as.matrix(
    raw_res[ , -1]
  )
  res_common = matrix( as.integer(res_common), ncol = dim(res_common)[2]  )
  
  ### 
  res_common = matrix( as.integer( unlist( sim_list ) ), ncol = length( sim_list)  )
  mapping    = match( rownames(res_common), vcf_fingerprint, nomatch = 0 )
  mapping[ mapping != 0 ] = 1
  
  adv = 0
  
  match_fp = function( sim_list_entry ){
    
    stat = round( (adv / as.double(nr_cls)) * 100, 1 )
    adv <<- adv + 1
    
    if ( stat != round( (adv / as.double(nr_cls)) * 100, 1 ) )
      print( paste( round( (adv / as.double(nr_cls)) * 100, 1 ), "% finished", sep =" " ) )
    
    return( sum( as.integer(mapping) & as.integer( sim_list_entry ) ) )
  }
  identification = lapply( sim_list, FUN = match_fp )
  
  ###
  
  hits            = colSums(res_common)
  
  match_index     = hits != 0
  candidates      = cl_names[ order(hits, decreasing = T)  ]
  
  print( paste0("Best candidate: ", candidates[1] )  )

  #res_common_filt = res_common[ match_index ]
  #table( res_common_filt)
  
  #identify_vcf_fingerprint( vcf_fingerprint  )
    
}