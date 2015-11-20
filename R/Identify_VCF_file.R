#' Parses the vcf file and predicts the identity of the sample
#' @export
identify_vcf_file = function( vcf_file_path ){
  
  print( paste0( "Creating fingerprint from VCF file ", vcf_file_path  ) )
  vcf_fingerprint = parse_vcf_file( vcf_file_path )
  
  sim_list_file = paste( system.file("", package="Younikorn"), "simlist.RData", sep = "/")
  print( paste0( "Loading similarity data from file ",  sim_list_file )  )
  
  attach( sim_list_file  )
  
  print( "Finished loading similarity data"  )
  
  #init_connection = function ( db_path ){

   # if ( grepl( "Younikorn.db", c(db_path)) != T )
    #  db_path = paste( db_path, "inst/Younikorn.db", sep = "/" )

  #  drv = dbDriver("SQLite")
  #  con = dbConnect(
   #   drv,
  #    dbname = db_path
  #  )
  #  con = src_sqlite( con@dbname )
      
  #  return( con )
  #}
  
  #full_con = init_connection( db_path )
  
  #raw_res = as.data.frame(

  #   dplyr::filter( 
  #    tbl( full_con, from = "similarity_matrix" ),
  #    mutational_similarity_marker %in% vcf_fingerprint
 #    )
 # )
 # cl_names = colnames(raw_res)[-1]

 # res_common = as.matrix(
 #   raw_res[ , -1]
 # )
 # res_common = matrix( as.integer(res_common), ncol = dim(res_common)[2]  )
  
  ###
  #res_common = matrix( as.integer( unlist( sim_list ) ), ncol = length( sim_list)  )
  #mapping    = match( rownames(res_common), vcf_fingerprint, nomatch = 0 )
  #mapping[ mapping != 0 ] = 1
  
  adv = 0
  nr_cls = length( sim_list)
  
  match_fp = function( sim_list_entry ){
    
    stat = round( (adv / as.double(nr_cls)) * 100, 1 )
    adv <<- adv + 1
    
    if ( stat != round( (adv / as.double(nr_cls)) * 100, 1 ) )
      print( paste( round( (adv / as.double(nr_cls)) * 100, 1 ), "% finished", sep =" " ) )
    
    return( sum( as.integer(mapping) & as.integer( sim_list_entry ) ) )
  }
  
  hits = unlist( lapply( sim_list, FUN = match_fp ) )
  hits = unlist(identification)
  
  match_index     = hits != 0
  candidates      = cl_data$CL[ order(hits, decreasing = T)  ]

  print( paste0("Best candidate: ", candidates[1] )  )

  #res_common_filt = res_common[ match_index ]
  #table( res_common_filt)
  
  #identify_vcf_fingerprint( vcf_fingerprint  )
    
}