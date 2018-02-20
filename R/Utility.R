write_mutation_grange_objects = function(
    g_mat,
    library_name,
    mutational_weight_inclusion_threshold,
    ref_gen,
    type = ""
){
    
    package_path = system.file("", package = "Uniquorn")
    
    rdata_path =  paste( c(
        package_path,"/Libraries/",
        ref_gen,"/",library_name,
        "/W",as.character(mutational_weight_inclusion_threshold),
        "_Uniquorn_DB.RData"),
        sep ="", collapse= ""
    )
    
    if ( type != ""){
        rdata_path = str_replace(rdata_path, pattern = "_Uniquorn_DB.RData",
            paste( c( ".",type,"_Uniquorn_DB.RData"), sep = "",collapse = "" )
        )
    }
    
    dir.create(
        paste( 
            c(package_path,"/Libraries/",ref_gen,"/"),
            collapse = "",
            sep =""
        ),
        showWarnings = FALSE)
    dir.create(
        paste( 
            c(package_path,"/Libraries/",ref_gen,"/",library_name),
            collapse = "",
            sep =""
        ),
        showWarnings = FALSE
    )
    saveRDS(g_mat,rdata_path)
}

#' read_mutation_grange_objects
#' 
#' Read the GRange object for a specific library
#' 
#' @param library_name a character string giving the name of the library
#' @param mutational_weight_inclusion_threshold a numerical giving the
#'  lower bound for mutational weight to be included
#' @param ref_gen Reference genome version. All training sets are 
#'  associated with a reference genome version. Default: GRCH37
#' @param test_mode Is this a test? Just for internal use
#' @importFrom IRanges IRanges
#' @return The R Table sim_list which contains the CoSMIC CLP fingerprints 
read_mutation_grange_objects = function(
    library_name,
    mutational_weight_inclusion_threshold,
    ref_gen,
    test_mode
){
    
    package_path = system.file("", package = "Uniquorn")
    
    rdata_path =  paste( c(
        package_path,"/Libraries/",
        ref_gen,"/",library_name,
        "/W",as.character(mutational_weight_inclusion_threshold),
        "_Uniquorn_DB.RData"),
        sep = "",
        collapse = ""
    )
    
    if (! file.exists( rdata_path )){
        message(paste(c("Database not found: ", rdata_path,
                        ", creating new database."),collapse = "",sep =""))
        g_mat = GenomicRanges::GRanges(
            seqnames = NULL,
            IRanges(
                start = NULL,
                end = NULL
            )
        )
        
        if (test_mode == FALSE){
            dir.create(
                paste( 
                    c(package_path,"/Libraries/"),
                    collapse = "",
                    sep = ""
                ),
                showWarnings = FALSE
            )
            dir.create(
                paste( 
                    c(package_path,"/Libraries/",ref_gen,"/"),
                    collapse = "",
                    sep = ""
                ),
                showWarnings = FALSE
            )
            dir.create(
                paste( 
                    c(package_path,"/Libraries/",ref_gen,"/",library_name),
                    collapse = "",
                    sep = ""
                ),
                showWarnings = FALSE
            )
        }
    } else {
        g_mat = readRDS(rdata_path)
    }
    
    return(g_mat)
}
