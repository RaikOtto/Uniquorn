#' show_contained_ccls
#' 
#' This function displays the names, amount of mutations and the overall 
#' weight of the mutations of all contained cancer cell line fingerprints 
#' for a chosen reference genome and optional library.
#' 
#' @param ref_gen a character vector specifying the reference genome version.
#'  All training sets are associated with a reference genome version.
#'  Default is \code{"GRCH37"}.
#' @param verbose Should DB informations be printed?
#' @return R table which contains identifiers of all cancer cell line samples
#'  which match the specified parameters (reference genome and library).
#' @usage 
#' show_contained_ccls(ref_gen, verbose)
#' @examples
#' ## Show all contained cancer cell lines for reference GRCH37:
#' show_contained_ccls(ref_gen = "GRCH37", verbose = TRUE)
#' @import GenomicRanges stringr
#' @export
show_contained_ccls = function(
    ref_gen = "GRCH37",
    verbose = TRUE
){
    
    package_path = system.file("", package = "Uniquorn")
    library_path =  paste( c( package_path,"/Libraries/",ref_gen),
        sep ="", collapse= "")
    
    if ( ! dir.exists( library_path ))
        stop("No libraries found!")
    
    libraries = list.dirs(library_path, full.names = FALSE)
    libraries = libraries[ libraries != ""]
    
    ccls_all <<- data.frame(
        "CCL" = as.character(),
        "Library" = as.character(),
        "W0" = as.character(),
        "W25" = as.character(),
        "W05" = as.character(),
        "W1" = as.character()
    )
    
    for (library_name in libraries){
        
        stats_path = paste( c( library_path,"/",library_name,
            "/CCL_List_Uniquorn_DB.RData"), sep ="", collapse= "")
        g_library = readRDS(stats_path)
        cl_id_no_library = str_replace_all( g_library$CCL, pattern = 
            paste("_", library_name, sep = ""), "" )
        g_library$CCL = cl_id_no_library
        
        if (verbose) {
            message(library_name," amount CCLs: ",
                as.character(nrow(g_library)))
        }
        g_library$Library = rep(library_name, nrow(g_library))
        g_library = g_library[c("CCL","Library","W0","W25","W05","W1")]
        ccls_all <<- rbind(ccls_all, g_library )
    }
    return(ccls_all)
}


#' All variants contained in reference library
#' 
#' This function shows all variants contained in a reference library
#' for a given inclusion weight. Default inclusion weight is 0 
#' (all variants).
#' 
#' @param ref_gen a character vector specifying the reference genome 
#' version. All training sets are associated with a reference genome version. 
#' Default is \code{"GRCH37"}.
#' @param library_name Name of the reference library.
#' @param mutational_weight_inclusion_threshold Include only mutations 
#' with a weight of at least x. Range: 0.0 to 1.0. 1 = unique to CL. 
#' ~0 = found in many CL samples. 
#' @usage
#' show_contained_variants_in_library(
#'     ref_gen,
#'     library_name,
#'     mutational_weight_inclusion_threshold
#' )
#' @import stringr
#' @examples
#' ## Show all variants contained in reference library CELLMINER
#' show_contained_variants_in_library(
#'     ref_gen = "GRCH37",
#'     library_name = "CELLMINER",
#'     mutational_weight_inclusion_threshold = 0
#' )
#' @return Returns a GenomicRanges object that contains the variants
#' @export
show_contained_variants_in_library = function(
    ref_gen = "GRCH37",
    library_name,
    mutational_weight_inclusion_threshold = 0
){
    
    message("Entered reference genome: ", ref_gen, 
        ". Entered library name: ", library_name)
    
    package_path = system.file("", package = "Uniquorn")
    library_path =  paste(
        c(package_path, "/Libraries/", ref_gen), sep = "", collapse= ""
    )
    
    if ( ! dir.exists( library_path ))
        stop("No libraries found!")
    
    libraries = list.dirs(library_path,full.names = FALSE)
    libraries = libraries[ libraries != ""]
    
    if (! (library_name %in% libraries ))
        stop("Could not find library!")
    
    g_mat = read_mutation_grange_objects(
        library_name = library_name,
        ref_gen = ref_gen,
        mutational_weight_inclusion_threshold = 
            mutational_weight_inclusion_threshold,
        test_mode = FALSE
    )
    
    return(g_mat)
}

#' Variants In Cancer Cell Line
#' 
#' This function shows all mutations present in the database 
#' for a selected cancer cell line and reference genome.
#' 
#' @param name_ccl a character vector giving the identifier 
#' of the cancer cell line for which mutations will be shown.
#' @param ref_gen a character vector specifying the reference 
#' genome version. All training sets are associated with a 
#' reference genome version. Default is \code{"GRCH37"}.
#' @param library_name Name of the reference library 
#' @param mutational_weight_inclusion_threshold Include only mutations 
#' with a weight of at least x. Range: 0.0 to 1.0. 1= unique to CL. 
#' ~0 = found in many CCL samples. 
#' @import stringr
#' @usage 
#' show_contained_variants_for_ccl(
#'     name_ccl, 
#'     ref_gen,
#'     library_name,
#'     mutational_weight_inclusion_threshold
#' )
#' @examples 
#' ## Show all mutations for Cancer Cell Line 'SK_OV_3'
#' show_contained_variants_for_ccl(
#'     name_ccl = "SK_OV_3",
#'     ref_gen = "GRCH37",
#'     library_name = "CELLMINER",
#'     mutational_weight_inclusion_threshold = 0
#' )
#' @return GenomicRanges object that contains the ccl's variants
#' @export
show_contained_variants_for_ccl = function(
    name_ccl,
    ref_gen = "GRCH37",
    library_name,
    mutational_weight_inclusion_threshold = 0
){
    
    message("Entered reference genome: ", ref_gen, 
        ". Entered library name: ", library_name,
        ". Entered name of ccl: ", name_ccl)
        
    package_path = system.file("", package = "Uniquorn")
    library_path =  paste(c(package_path, "/Libraries/", ref_gen),
                            sep ="", collapse= "")
    
    if ( ! dir.exists( library_path ))
        stop("No libraries found!")
    
    libraries = list.dirs(library_path,full.names = FALSE)
    libraries = libraries[ libraries != ""]
    
    if (! (library_name %in% libraries ))
        stop("Could not find library!")
    
    g_mat = read_mutation_grange_objects(
        library_name = library_name,
        ref_gen = ref_gen,
        mutational_weight_inclusion_threshold = 
            mutational_weight_inclusion_threshold,
        test_mode = FALSE
    )
    
    g_query_index = stringr::str_detect(g_mat$Member_CCLs, pattern = name_ccl)
    g_query = g_mat[g_query_index,]
    
    return(g_query)
}

#' Cancer cell lines with specific variant
#' 
#' This function displays all cancer cell lines in the database which 
#' contain a specified variant. Utilizes closed interval coordinates.
#' 
#' @param start Start coordinate
#' @param end Stop coordinate
#' @param chromosome Chromosome, 'chr' prefixes are ignored
#' @param ref_gen a character vector specifying the reference genome 
#' version. All training sets are associated with a reference genome 
#' version. Default is \code{"GRCH37"}.
#' @param library_name Name of the reference library 
#' @param mutational_weight_inclusion_threshold Include only mutations 
#' with a weight of at least x. Range: 0.0 to 1.0. 1= unique to CL. 
#' ~0 = found in many CCL samples. 
#' @import stringr
#' @importFrom IRanges IRanges
#' @importFrom IRanges subsetByOverlaps
#' @usage 
#' show_which_ccls_contain_variant(
#'     start,
#'     end,
#'     chromosome,
#'     ref_gen,
#'     library_name,
#'     mutational_weight_inclusion_threshold
#' )
#' @examples
#' show_which_ccls_contain_variant(
#'     start = 92030762,
#'     end = 92030762,
#'     chromosome = 8,
#'     ref_gen = "GRCH37",
#'     library_name = "CELLMINER",
#'     mutational_weight_inclusion_threshold = 0
#' )
#' @return Returns a GenomicRanges object that contains the variant if present.
#' Member ccls can be found in the $Member_ccl vector
#' @export
show_which_ccls_contain_variant = function(
    start,
    end,
    chromosome,
    ref_gen = "GRCH37",
    library_name,
    mutational_weight_inclusion_threshold = 0
){
    
    message("Entered reference genome: ", ref_gen, 
        ". Entered library name: ", library_name)
    
    package_path = system.file("", package = "Uniquorn")
    library_path = paste(
        c(package_path, "/Libraries/", ref_gen), sep ="", collapse= ""
    )
    
    if ( ! dir.exists( library_path ))
        stop("No libraries found!")
    
    libraries = list.dirs(library_path,full.names = FALSE)
    libraries = libraries[ libraries != ""]
    
    if (! (library_name %in% libraries ))
        stop("Could not find library!")
    
    g_mat = read_mutation_grange_objects(
        library_name = library_name,
        ref_gen = ref_gen,
        mutational_weight_inclusion_threshold = 
        mutational_weight_inclusion_threshold,
        test_mode = FALSE
    )
    
    chromosome = stringr::str_to_upper(chromosome)
    chroms = stringr::str_replace( chromosome, pattern = "^CHR" , "")
    
    g_query = GenomicRanges::GRanges(
        seqnames = c( chroms ),
        IRanges::IRanges(
            start = as.integer( c(start) ),
            end = as.integer( c(end) )
        )
    )
    
    g_query = IRanges::subsetByOverlaps(g_mat, g_query)
    return(g_query)
}

#' Library Name Reader
#' 
#' This function procides information on the reference library names
#' 
#' @param ref_gen a character vector specifying the reference genome 
#' version. All training sets are associated with a reference genome 
#' version. Default is \code{"GRCH37"}.
#' @import stringr
#' @usage
#' read_library_names(ref_gen)
#' @examples
#' read_library_names(ref_gen = "GRCH37")
#' @return Returns a character vector of the contained libraries
#' @export
read_library_names = function(
    ref_gen
){
    
    package_path = system.file("", package = "Uniquorn")
    library_path =  paste(c(package_path, "/Libraries/", ref_gen, "/"), 
        sep ="", collapse= "")
    
    library_names = list.dirs(library_path, full.names = FALSE)
    library_names = library_names[library_names!= ""]
    
    return(library_names)
}
