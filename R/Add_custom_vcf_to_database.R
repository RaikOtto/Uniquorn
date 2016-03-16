#' Adds a custom vcf file to the three existing cancer cell line panels
#' @param vcf_file_path Input vcf file. Only one sample column allowed.
#' @param ref_gen Reference genome version. All training sets are associated 
#' with a reference genome version. Default: GRCH37
#' @param name_cl Name of the to-be-added cancer cell line sample. '_CUSTOM' 
#' will automatically be added as suffix.
#' @param safe_mode Only add mutations to the database where there already are 
#' mutations found in the cannonical cancer cell lines. This is a safety 
#' mechanism against overfitting if there are too few custom training samples.
#' @param distinct_mode Show training data for the commonly or separately 
#' normalized training sets. Options: TRUE/ FALSE
#' @param test_mode Is this a test? Just for internal use
#' @return Message if the adding has succeeded
#' @import DBI
#' @usage 
#' add_custom_vcf_to_database( 
#'  vcf_file_path, 
#' ref_gen = "GRCH37", 
#' name_cl = "", 
#' safe_mode = FALSE, 
#' distinct_mode = TRUE,
#' test_mode = FALSE)
#' @examples 
#' HT29_vcf_file = system.file("extdata/HT29.vcf.gz", package="Uniquorn");
#' add_custom_vcf_to_database( 
#' vcf_file_path = HT29_vcf_file,
#' name_cl = "",
#' ref_gen = "GRCH37",
#' safe_mode = FALSE,
#' distinct_mode = TRUE,
#' test_mode = TRUE )
#' @export
add_custom_vcf_to_database = function( 
    vcf_file_path,
    ref_gen = "GRCH37",
    name_cl = "",
    safe_mode = FALSE,
    distinct_mode = TRUE,
    test_mode = FALSE
    ){
    
    # pre processing

    name_cl = stringr::str_to_upper(name_cl)
    
    base::print( base::paste0( "Reference genome: ",ref_gen ) )
    
    sim_list_stats = initiate_db_and_load_data( ref_gen = ref_gen, 
        distinct_mode = distinct_mode, request_table = "sim_list_stats")
    
    sim_list_stats = sim_list_stats[ sim_list_stats$Ref_Gen == ref_gen,]
    
    if ( dim(sim_list_stats)[1] == 0 )
        warning( paste( "Warning: Identification might be spurious due to 
        low amount of training sample. No cancer cell line data stored 
        at this point for reference genome:", ref_gen, sep =" " ) )
    
    # reading vcf
    
    if ( name_cl == "" ){
        
        vcf_identifier = as.character( 
            utils::tail( 
                unlist( 
                    stringr::str_split( 
                        vcf_file_path, 
                        "/" 
                    )
                ), 
            1)
        )
        name_cl = utils::head( 
            unlist( stringr::str_split( vcf_identifier, ".vcf|.VCF" ) ) ,
        1)
        name_cl = stringr::str_to_upper( base::paste( 
            name_cl, "CUSTOM", sep ="_"  
        ))
        print( base::paste0( 
            "No cl name provided, adding auto-generated fingerprint: ",
            name_cl
        ) )
        
    } else {
        
        name_cl = stringr::str_to_upper( base::paste(name_cl, "custom",
            sep ="_"
        ) )
        base::print( base::paste0( 
            "Adding fingerprint with user-defined name: ",
        name_cl ) )
    }
    
    name_present = base::grepl( name_cl, sim_list_stats$CL )
    
    if( sum( name_present ) > 0 ){ 
        stop( paste0( c("Fingerprint with name ",name_cl, 
            " already present in database. Please change the 
            name or remove the old cancer cell line."), 
            collapse = "" )  )
    }
    
    sim_list = initiate_db_and_load_data( ref_gen = ref_gen, 
        distinct_mode = distinct_mode, request_table = "sim_list" )
    sim_list = sim_list[ sim_list$Ref_Gen == ref_gen,]
    sim_list = sim_list[, base::which( 
        colnames(sim_list) != "Ref_Gen"
    ) ]
    sim_list = sim_list[, base::which( 
        colnames(sim_list) != "Weight"  )]
    
    print( paste0( "Building fingerprint from file ",
        vcf_file_path )  )
    
    vcf_fingerprint = as.character( parse_vcf_file(
        vcf_file_path ) )
    
    if( safe_mode )
        vcf_fingerprint = vcf_fingerprint[ which( 
        vcf_fingerprint %in% sim_list$Fingerprint ) ]
    
    vcf_fingerprint = data.frame( 
        "Fingerprint" = vcf_fingerprint,
        "CL"  = rep( name_cl, length(as.character(vcf_fingerprint)) )
    )
    
    sim_list = rbind( sim_list, vcf_fingerprint)
    
    print("Finished parsing, aggregating over parsed Cancer Cell Line data")
    
    res_vec = re_calculate_cl_weights( sim_list = sim_list, 
        ref_gen = ref_gen, distinct_mode = TRUE )
    
    print("Finished aggregating, saving to database")
    
    write_data_to_db( content_table = res_vec[1], "sim_list",       
        ref_gen = "GRCH37",  distinct_mode = distinct_mode, 
        overwrite = TRUE, test_mode = test_mode )
    
    write_data_to_db( content_table = res_vec[2], "sim_list_stats", 
        ref_gen = "GRCH37", distinct_mode = distinct_mode, 
        overwrite = TRUE, test_mode = test_mode )
    
    print("Finished adding CL")
}