#' Filter Parsed VCF Files
#' 
#' Intern utility function. Filters the parsed VCF file for all informations
#' except for the start and length of variations/mutations.
#' 
#' @param vcf_file character string giving the path to the vcf file
#'  on the operating system.
#' @param ref_gen Reference genome version
#' @param library_name Name of the reference library
#' @usage 
#' parse_vcf_file(
#'     vcf_file,
#'     ref_gen,
#'     library_name
#' )
#' @import GenomicRanges stringr
#' @importFrom IRanges IRanges
#' @importFrom VariantAnnotation readVcf
#' @return Loci-based DNA-mutational fingerprint of the cancer cell line
#'  as found in the input VCF file.
parse_vcf_file = function(
    vcf_file,
    ref_gen,
    library_name
){
    switch(ref_gen,
        "GRCH37" = {ref_gen = "hg19"}
    )
    
    g_query =  VariantAnnotation::readVcf(file = vcf_file)
    
    # process variants
    chroms = stringr::str_replace(as.character(seqnames(g_query)),
                         pattern = "chr|CHR", "")
    start_var = start(g_query)
    end_var = start_var + width(g_query) - 1
    
    chroms_pure = grep(chroms, pattern = "_", invert = TRUE)
    chroms      = chroms[chroms_pure]
    chroms      = str_replace(chroms, pattern = "^chr", "")
    start_var   = start_var[chroms_pure]
    end_var     = end_var[chroms_pure]
    
    cl_id = gsub("^.*/", "", vcf_file)
    cl_id = gsub(".vcf", "", cl_id, fixed = TRUE)
    cl_id = gsub(".VCF", "", cl_id, fixed = TRUE)
    cl_id = gsub(".hg19", "", cl_id, fixed = TRUE)
    cl_id = toupper(cl_id)
    cl_id = str_replace_all(cl_id, pattern = "\\.", "_")
    cl_id = paste0(cl_id, "_", library_name)
    
    # build fingerprint and return
    g_query = GenomicRanges::GRanges(
        seqnames = chroms,
        ranges = IRanges(start = start_var, end = end_var),
        Member_CCLs = rep(cl_id, length(g_query))
    )

    return(g_query)
}

write_w0_and_split_w0_into_lower_weights = function(
    g_mat,
    library_name,
    ref_gen
){
    
    ccl_stats <<- data.frame(
        "CCL" = names(table(as.character(
            unlist(str_split(g_mat$Member_CCLs, pattern = ","))))),
        "Mut_count_0" = as.character(table(as.character(unlist(
            str_split(g_mat$Member_CCLs, pattern = ",")))))
    )
    
    if ( "Member_CCLs" %in% names(GenomicRanges::mcols(g_mat)) ){
        
        member_count = as.integer( str_count(g_mat$Member_CCLs, pattern = ","))
        member_keep_index    = which( member_count < 10 )
        member_exclude_index = which( member_count >= 10 )
        
        g_mat_exclude <<- g_mat[member_exclude_index]
        g_mat         <<- g_mat[member_keep_index]
        
        write_mutation_grange_objects(
            mutational_weight_inclusion_threshold = 0,
            g_mat = g_mat_exclude,
            library_name = library_name,
            ref_gen = ref_gen,
            type = "excluded_variants"
        )
    }
    
    write_mutation_grange_objects(
        mutational_weight_inclusion_threshold = 0,
        g_mat = g_mat,
        library_name = library_name,
        ref_gen = ref_gen
    )
    
    for(mwit in c(.25,.5,1.0)){
        
        if ( !(  "Member_CCLs" %in% names(GenomicRanges::mcols(g_mat)) )){
            
            hit_index = 0
            
        } else {
            
            hit_index = which( str_count(
                GenomicRanges::mcols(g_mat)$Member_CCLs, pattern = ","
                ) <= as.integer( round(1/mwit) ) )
        }
        
        mwit_g_mat = g_mat[hit_index]
        
        write_mutation_grange_objects(
            mutational_weight_inclusion_threshold = mwit,
            g_mat = mwit_g_mat,
            library_name = library_name,
            ref_gen = ref_gen
        )
        
        mut_t = table(as.character(unlist(str_split( 
            GenomicRanges::mcols(mwit_g_mat)$Member_CCLs,pattern = "," ))))
        
        ccl_count = rep("0", nrow(ccl_stats))
        ccl_match = match( ccl_stats$CCL, names(mut_t), nomatch = 0 )
        ccl_count[ccl_match] = mut_t[ccl_match != 0]
        ccl_stats = cbind(ccl_stats, ccl_count)
    }
    colnames(ccl_stats) = c("CCL","W0","W25","W05","W1")
    
    package_path = system.file("", package = "Uniquorn")
    
    input_ccls = unique(as.character(unlist(
        str_split( g_mat$Member_CCLs, pattern = "," ) ) ))
    output_ccls = unique(as.character(unlist(
        str_split( ccl_stats$CCL, pattern = "," ) ) ))
    missing_ccls = sum((input_ccls %in%  output_ccls) == FALSE)
    
    if( missing_ccls != 0) stop("Loading of CCLs has not worked!")
    
    CCL_stats_data_path =  paste( c(
        package_path,"/Libraries/",
        ref_gen,"/",library_name,
        "/CCL_List_Uniquorn_DB.RData"),
        sep ="", collapse= ""
    )
    saveRDS(ccl_stats,CCL_stats_data_path)
}