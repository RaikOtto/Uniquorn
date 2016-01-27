#' initiate_canonical_databases
#' 
#' Parses data into r list variable
#' 
#' @param cosmic_file The path to the cosmic DNA genotype data file. Ensure that the right reference genome is used
#' @param ccle_file The path to the ccle DNA genotype data file. Ensure that the right reference genome is used
#' @param ref_gen Reference genome version
#' @param distinct_mode Should the mutational weights be calculated for all panels together or each for itelf? Recommendation: Seperately
#' @import DBI R.utils RSQLite
#' @usage 
#' initiate_canonical_databases( 
#' cosmic_file = "CosmicCLP_MutantExport.tsv", 
#' ccle_file = "CCLE_hybrid_capture1650_hg19_NoCommonSNPs_CDS_2012.05.07.maf", 
#' ref_gen = "GRCH37",
#' distinct_mode = TRUE)
#' @export
initiate_canonical_databases = function(
    cosmic_file = "CosmicCLP_MutantExport.tsv",
    ccle_file = "CCLE_hybrid_capture1650_hg19_NoCommonSNPs_CDS_2012.05.07.maf",
    ref_gen = "GRCH37",
    distinct_mode = TRUE
    ){

    print( c( "Reference genome: ", ref_gen )  )
    
    ### pre-processing
    
    package_path  = system.file("", package="Uniquorn")
    database_path = base::paste( package_path, "uniquorn_distinct_panels_db.sqlite", sep ="/" )
    database_default_path = base::paste( package_path, "uniquorn_db_default.sqlite", sep ="/" )
    
    if (!distinct_mode)
        database_path   =  base::paste( package_path, "uniquorn_non_distinct_panels_db.sqlite", sep ="/" )

    if ( base::file.exists(database_path) )
        base::file.copy( from = database_path, to = database_default_path, overwrite = TRUE )
    
    sim_list       = inititate_db_and_load_data( ref_gen = ref_gen, distinct_mode = distinct_mode, request_table = "sim_list" )
    
    sim_list = sim_list[, which( colnames(sim_list) != "Ref_Gen"  ) ]
    sim_list = sim_list[, which( colnames(sim_list) != "Weight"  ) ]

    parse_files = c()
  
    if (file.exists(cosmic_file)){
      
        print( c( "Found CoSMIC: ", file.exists(cosmic_file) )  )
      
        if ( grepl( ".gz$", stringr::str_to_lower( cosmic_file ) ) ){

            gunzip( cosmic_file, overwrite = TRUE )
        }
        cosmic_file = stringr::str_replace( cosmic_file, ".gz$|.GZ$", "" )
      
        sim_list = parse_cosmic_genotype_data( cosmic_file, sim_list )
    }
  
    if (file.exists(ccle_file)){
      
      print( c( "Found CCLE: ", file.exists( ccle_file ) )  )
      sim_list = parse_ccle_genotype_data( ccle_file, sim_list )
    }
    
    if (length(parse_files) == 0)
      stop("Did not find CCLE & CoSMIC CLP file! Aborting.")
    
    # overwrite existing db
    if (file.exists(database_path))
        file.remove( database_path )
    
    print("Finished parsing, aggregating over parsed Cancer Cell Line data")

    res_vec = re_calculate_cl_weights( sim_list = sim_list, ref_gen = ref_gen, distinct_mode = TRUE )
  
    print("Finished aggregating, saving to database")
    
    write_data_to_db( content_table = res_vec[1], "sim_list",       ref_gen = "GRCH37", distinct_mode = distinct_mode, overwrite = TRUE )
    write_data_to_db( content_table = res_vec[2], "sim_list_stats", ref_gen = "GRCH37", distinct_mode = distinct_mode, overwrite = TRUE )
    
    print ("Initialization of Uniquorn DB finished")
}