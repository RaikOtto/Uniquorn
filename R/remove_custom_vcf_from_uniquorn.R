#' Removes a cancer cell line training fingerprint (vcf file) from the database. The names of all training sets can 
#' be seen by using the function \code{show_contained_cls}.
#' @param name_cl name of the cancer cell line training fingerprintt
#' @param ref_gen Reference genome version. All training sets are associated with a reference genome version. Default: GRCH37
#' @param test_mode Is this a test? Just for internal use
#' @import DBI
#' @usage 
#' remove_custom_vcf_from_database( 
#' name_cl, 
#' ref_gen = "GRCH37", 
#' test_mode = FALSE)
#' @examples 
#' remove_custom_vcf_from_database( 
#' name_cl = "HT29_CELLMINER", 
#' ref_gen = "GRCH37",
#' test_mode = TRUE )
#' @return Message that indicates if the removal was succesful
#' @export
remove_custom_vcf_from_database = function( 
    name_cl,
    ref_gen = "GRCH37",
    test_mode = FALSE
){
    
    # pre processing

    name_cl = stringr::str_to_upper(name_cl)
    
    print(paste0("Reference genome: ",ref_gen))
    
    sim_list       = initiate_db_and_load_data( ref_gen = ref_gen, request_table = "sim_list" )
    sim_list_stats = initiate_db_and_load_data( ref_gen = ref_gen, request_table = "sim_list_stats" )

    if ( sum( grepl( name_cl, sim_list_stats$CL ) ) == 0 ){
        
        stop(paste0("No training set for a cancer cell line found for the name: ", name_cl))
        
    } else {
        
        print( 
            paste0(
                c("Found CL ",name_cl,
                ". Removing from database and recalculating training-sets."
            ), 
            collapse = "") 
        )
               
    }
    
    sim_list = sim_list[ sim_list$CL != name_cl,] # exclude sample here
    sim_list = sim_list[, which( colnames(sim_list) != "Ref_Gen"  ) ]
    sim_list = sim_list[, which( colnames(sim_list) != "Weight"  ) ]
    
    print("Found & removed the sample. Re-calculating the Cancer cell line data")
    
    res_vec = re_calculate_cl_weights( sim_list = sim_list, ref_gen = ref_gen)
    
    print("Finished aggregating, saving to database")
    
    write_data_to_db( content_table = as.data.frame( res_vec[1] ), "sim_list",       ref_gen = "GRCH37", overwrite = TRUE, test_mode = test_mode )
    write_data_to_db( content_table = as.data.frame( res_vec[2] ), "sim_list_stats", ref_gen = "GRCH37", overwrite = TRUE, test_mode = test_mode )
    
    print(paste0( "Finished removing CL ", name_cl ))
}