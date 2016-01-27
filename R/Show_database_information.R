#' show_contained_cls
#' 
#' Show all cancer cell line identifier present in the database for a selected reference genome:
#' This function shows the names, amount of mutations/ variations, overall weight of the mutations of all contained training CLs 
#' for a chosen reference genome.
#' 
#' @param ref_gen Reference genome version. All training sets are associated with a reference genome version. Default: GRCH37
#' @param distinct_mode Show training data for the commonly or separately normalized training sets. Options: TRUE/ FALSE
#' @return R table which contains the identifier of all cancer cell line samples with the specific reference genome and the weight of all mutations
#' @usage 
#' show_contained_cls( 
#' ref_gen, 
#' distinct_mode )
#' @examples 
#' contained_cls = show_contained_cls( 
#' ref_gen = "GRCH37", 
#' distinct_mode = TRUE )
#' @import DBI RSQLite
#' @export
show_contained_cls = function( ref_gen = "GRCH37", distinct_mode = TRUE ){

    print(paste0("Reference genome: ",ref_gen))
    
    sim_list_stats = initiate_db_and_load_data( ref_gen = ref_gen, distinct_mode = distinct_mode, request_table = "sim_list_stats" )
    
    print( paste0( c("Found ", dim(sim_list_stats)[1], " many cancer cell lines fingerprints for reference genome ", ref_gen ), collapse = ""  )  )

    print( paste( "CoSMIC CLP: ", as.character( sum( grepl( "_COSMIC", sim_list_stats$CL ) ) ) ) )
    print( paste( "CCLE: ", as.character( sum( grepl( "_CCLE", sim_list_stats$CL ) ) ) ) )
    print( paste( "CellMiner: ", as.character( sum( grepl( "_CELLMINER", sim_list_stats$CL ) ) ) ) )
    print( paste( "CUSTOM: ", as.character( sum( grepl( "_CUSTOM", sim_list_stats$CL ) ) ) ) )
    
    return( sim_list_stats )
}

#' show_contained_mutations
#' 
#' Show all mutations present in the database for a selected reference Genome: 
#' This function shows all training-set mutations for a selected reference genome, i.e. the mutations that are being used
#' for identification of query cancer cell lines.
#' 
#' @param ref_gen Reference genome version
#' @param distinct_mode Show mutations for either distinct or non-distinct normalization of mutational weights
#' @usage 
#' show_contained_mutations( 
#' ref_gen, 
#' distinct_mode )
#' @examples 
#' contained_cls = show_contained_mutations( ref_gen = "GRCH37", distinct_mode = TRUE )
#' @return R Table which contains all mutations associated with a particular cancer cell line for a specified reference genome
#' @export
show_contained_mutations = function( ref_gen = "GRCH37", distinct_mode = TRUE ){
  
    print(paste0("Reference genome: ",ref_gen))
    
    sim_list = initiate_db_and_load_data( ref_gen = ref_gen, distinct_mode = distinct_mode, request_table = "sim_list" )
    
    print( paste0( c("Found ", dim(sim_list)[1], " many cancer cell lines associated mutations for reference genome ", ref_gen ), collapse = ""  )  )
  
    print( summary( sim_list ) )
  
    return( sim_list )  
}

#' show_contained_mutations_for_cl
#' 
#' Show all mutations present in the database for a selected cancer cell line and reference Genome
#' 
#' @param ref_gen Reference genome version
#' @param distinct_mode Show mutations for either distinct or non-distinct normalization of mutational weights
#' @param name_cl Name of the cancer cell line sample stored in the database
#' @import DBI
#' @usage 
#' show_contained_mutations_for_cl( 
#' name_cl, 
#' ref_gen, 
#' distinct_mode )
#' @examples 
#' SK_OV_3_CELLMINER_mutations = show_contained_mutations_for_cl(
#' name_cl = "SK_OV_3_CELLMINER_mutations",
#' ref_gen = "GRCH37",
#' distinct_mode = TRUE)
#' @return R table which contains all mutations associated with the defined cancer cell line and reference genome
#' @export
show_contained_mutations_for_cl = function( name_cl, ref_gen = "GRCH37", distinct_mode = TRUE){

    print(paste0("Reference genome: ",ref_gen))
    
    sim_list = initiate_db_and_load_data( ref_gen = ref_gen, distinct_mode = distinct_mode, request_table = "sim_list" )
  
    sim_list = sim_list[ sim_list$Ref_Gen == ref_gen,  ]
    mapping  = which( sim_list$CL %in% name_cl, arr.ind = TRUE  )
    sim_list = sim_list[ mapping,  ]
    
    if ( length( mapping ) == 0  ){
    
        message(paste0("Could not find the cancer cell line ",name_cl, " in the database."), collapse= "")
    
    } else {
    
        print( paste0( c("Found ", dim(sim_list)[1], " many mutations for cancer cell line", name_cl  ," for reference genome ", ref_gen ), collapse = ""  )  )
    }
    
    return(sim_list)
}


#' show_which_cls_contain_mutation
#' 
#' Show all cancer cell lines in the database which contained the specified mutation and reference Genome. Closed interval coordinates. Format mutation: CHR_START_STOP, e.g. 1_123_123
#' 
#' @param mutation_name Name of the mutation in the format CHROMOSOME_START_STOP, e.g. '11_244501_244510'
#' @param ref_gen Reference genome version
#' @param distinct_mode Show mutations for either distinct or non-distinct normalization of mutational weights
#' @usage 
#' show_which_cls_contain_mutation( 
#' mutation_name, 
#' ref_gen, 
#' distinct_mode)
#' @examples 
#' Cls_containing_mutations = show_which_cls_contain_mutation( 
#' mutation_name = "10_103354427_103354427", 
#' ref_gen = "GRCH37", 
#' distinct_mode = TRUE )
#' @import DBI
#' @return R table which contains all cancer cell line samples which contain the specified mutation with respect to the specified reference genome version
#' @export
show_which_cls_contain_mutation = function( mutation_name, ref_gen = "GRCH37", distinct_mode = TRUE){
  
    print(paste0("Reference genome: ",ref_gen))
    
    sim_list = initiate_db_and_load_data( ref_gen = ref_gen, distinct_mode = distinct_mode, request_table = "sim_list" )
  
    sim_list = sim_list[ sim_list$Ref_Gen == ref_gen,  ]
    mapping  = which( sim_list$Fingerprint %in% mutation_name, arr.ind = TRUE)
    sim_list = sim_list[ mapping,  ]
    
    if ( length( mapping ) == 0  ){
    
        message(paste0("Could not find any cancer cell line for the mutation ",mutation_name, " in the database."), collapse= "")
    
    } else {
    
        print( paste0( c("Found ", dim( sim_list )[1], " many cancer cell lines for mutation ", mutation_name  ," for reference genome ", ref_gen ), collapse = ""  )  )
    }
    
    return( sim_list )  
}