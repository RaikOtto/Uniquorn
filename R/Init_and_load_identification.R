#' init_and_load_identification
#' 
#' Initiate the analysis
#' Output basic information
#' 
#' \code{init_and_load_identification} parses vcf file and output 
#' basic information
#' @param verbose Print additional information
#' @param ref_gen Reference genome version. All training sets are 
#' associated with a reference genome version. Default: GRCH37
#' @param vcf_file Path to vcf_file
#' @param output_dir Output directory for identification results
#' @usage 
#' init_and_load_identification( 
#'     verbose,
#'     ref_gen,
#'     vcf_file,
#'     output_dir
#' )
#' @return Three file path instances and the fingerprint
init_and_load_identification = function(
    verbose,
    ref_gen,
    vcf_file,
    output_dir
){
    
    vcf_fingerprint = parse_vcf_file(vcf_file)
    vcf_file_name = tail(as.character(unlist(strsplit(vcf_file, "/"))), 1)
    
    if (output_dir == ""){
        output_file = base::paste(vcf_file, "uniquorn_ident.tab", sep = "_")
        
    } else if (base::dir.exists(output_dir)){
        output_file = base::paste(
            output_dir,
            base::paste( 
                vcf_file_name,
                "uniquorn_ident.tab",
                sep = "_"
            ),
            sep = "/"
        ) 
    }
    output_file_xls = base::gsub(".tab$", ".xls", output_file)
    
    res_list = list( 
        "output_file"     = output_file,
        "output_file_xls" = output_file_xls,
        "vcf_fingerprint" = vcf_fingerprint,
        "vcf_file_name"   = vcf_file_name
    )
    
    return(res_list)
}