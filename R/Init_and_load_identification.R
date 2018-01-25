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
#' @param output_file Path to output report file
#' @param n_threads Specifies number of threads to be used
#' @import BiocParallel
#' @usage 
#' init_and_load_identification( 
#' verbose,
#' ref_gen,
#' vcf_file,
#' output_file,
#' n_threads)
#' @return Three file path instances and the fingerprint
init_and_load_identification = function( verbose, ref_gen, vcf_file, output_file, n_threads ){

    if ( verbose )  
        print( paste0("Assuming reference genome ", ref_gen) )
    
    ### pre processing
    
    if (verbose)
        print( paste0("Reading VCF file: ", vcf_file ) )
    
    vcf_fingerprint = parse_vcf_file( vcf_file, n_threads )
    vcf_file_name = utils::tail( as.character( unlist( stringr::str_split( vcf_file, "/" ) ) ), 1 )
    
    if ( output_file == ""  ){
        
        output_file = base::paste( vcf_file, "uniquorn_ident.tab", sep ="_")
        
    }else if ( base::dir.exists( output_file ) ){
        
        output_file = base::paste(
            output_file,
            base::paste( 
                vcf_file_name,
                "uniquorn_ident.tab",
                sep ="_"
            ),
            sep = "/"
        ) 
    }
    output_file_xls = stringr::str_replace( output_file, ".tab$", ".xls" ) 
    
    if ( verbose )
        print( "Finished reading the VCF file, loading database" )

    res_list = list( 
        "output_file"     = output_file,
        "output_file_xls" = output_file_xls,
        "vcf_fingerprint" = vcf_fingerprint,
        "vcf_file_name"   = vcf_file_name
    )
    
    return( res_list )
}