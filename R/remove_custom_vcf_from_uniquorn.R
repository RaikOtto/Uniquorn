#' Remove Cancer Cell Line
#'
#' This function removes a cancer cell line training fingerprint (VCF file)
#'  from the database. The names of all training sets can 
#'  be seen by using the function \code{show_contained_cls}.
#' 
#' @param ccl_names A character vector giving the names of the cancer cell line
#'  identifiers to be removed. Can be one or many
#' @param ref_gen A character vector specifying the reference genome version.
#'  All training sets are associated with a reference genome version. 
#'  Default is \code{"GRCH37"}.
#' @param library_name Name of the library from which ccls are to be removed
#' @param test_mode Signifies if this is a test run
#' @import GenomicRanges stringr
#' @usage 
#' remove_ccls_from_database(ccl_names, ref_gen = "GRCH37",
#'     library_name, test_mode = FALSE)
#' @examples 
#' remove_ccls_from_database(
#'     ccl_names = "HT29",
#'     ref_gen = "GRCH37",
#'     library_name = "CUSTOM",
#'     test_mode = TRUE
#' )
#' @return Message that indicates whether the removal was succesful.
#' @export
remove_ccls_from_database = function( 
    ccl_names,
    ref_gen = "GRCH37",
    library_name,
    test_mode = FALSE
){
    message("Reference genome: ", ref_gen)
    
    package_path = system.file("", package = "Uniquorn")
    library_path = paste(c(package_path, "/Libraries/", ref_gen),
                            sep ="", collapse= "")
    
    if (length(library_name) > 1)
        stop("Cannot process more than one library per run, please 
            provide single libraries.")
    
    g_mat = read_mutation_grange_objects(
        library_name = library_name,
        ref_gen = ref_gen,
        mutational_weight_inclusion_threshold = 0,
        test_mode = test_mode
    )
    
    for (ccl_id in ccl_names){
        
        ccl_lib = paste(ccl_id, library_name, sep = "_")
        message("Deleting ", ccl_id)
        search_term = paste(c(
            paste( "(",ccl_lib,",",")", sep = "" ),
            paste( "(",ccl_lib,")", sep = "" ),
            "(.*,$)"
        ), collapse= "|")
        
        member_ccls = str_replace_all(mcols(g_mat)$Member_CCLs,
                pattern = search_term, replacement = "")
        member_ccls = str_replace_all(member_ccls, pattern = ",$", "")
        
        non_empty_vec = member_ccls != ""
        g_mat = g_mat[non_empty_vec,]
        member_ccls = member_ccls[ non_empty_vec ]
        mcols(g_mat)$Member_CCLs = as.character(member_ccls)
        
        message(
            "Excluded ",
            as.character(sum(non_empty_vec == FALSE)),
            " variant entries after removal."
        )
        
        if( test_mode == FALSE ){
            stats_path = paste( c( library_path,"/",library_name,
                "/CCL_List_Uniquorn_DB.RData"), sep ="", collapse= "")
            g_library = readRDS(stats_path)
            g_library = data.frame(g_library, stringsAsFactors = FALSE)
            g_library = g_library[ g_library$CCL!= "",]
            g_library = g_library[!is.na(g_library$CCL),]
            g_library = g_library[g_library$CCL != ccl_id,]
            saveRDS(g_library, file = stats_path)
        }
    }
    
    #member_ccls = str_replace(member_ccls, pattern = ",$","")
    message("Finished removing all ccls. Recalculating DB")
    
    if( test_mode == FALSE){
        write_w0_and_split_w0_into_lower_weights(
            g_mat = g_mat,
            ref_gen = ref_gen,
            library_name = library_name
        )
    }
    message("Finished removing all cancer cell lines")
}


#' Remove entire Library from Database
#'
#' This function removes a entire library from the database by removing all 
#' associated cancer cell line fingerprints from the database.
#' 
#' @param library a character vector giving the names of the library to be
#'  removed.
#' @param ref_gen a character vector specifying the reference genome version.
#'  All training sets are associated with a reference genome version. 
#'  Default is \code{"GRCH37"}.
#' @param test_mode is this a test? Just for internal use.
#' @import stringr
#' @return Message that indicates whether the removal was succesful.
remove_library_from_database = function( 
    library,
    ref_gen = "GRCH37",
    test_mode = FALSE
){
    message("Reference genome: ", ref_gen)
    message("Removing library ", library, " and all associated",
        " cancer cell lines done.")
}