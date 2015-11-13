### parse files

#' Parses data into the sql database
#' Requires the currently (Nov 2015) github located package "chapmandu2/CancerCellLines"
#' @export
parse_data_into_younikorn_db = function( parser_path, db_path = system.file("database", package="Younikorn")  ){

  db_path  = paste0(db_path,"/Younikorn.db")
  print( paste0( "Parsing data and storing in db: ",db_path) )
  full_con = setupSQLite(db_path)

  # ids data
  
  idspath       = system.file("extdata", "CellLineIDNormalisationOct15.txt", package = "CancerCellLines")
  if (file.exists( idspath ) ){
    
    importCellLineIDs( idspath, full_con )
    print( paste0( "Parsed file ", idspath ) )
  }

  #infopath
  
  infopath = paste( parser_path, 'CCLE_sample_info_file_2012-10-18.txt', sep = "/" )
  
  if ( file.exists( infopath ) ){
    
    importCCLE_info(infopath , full_con)
    print(paste0("Parsed file ", infopath))
  }
  
  # ccle genotype data

  hybcappath = paste( parser_path, 'CCLE_hybrid_capture1650_hg19_NoCommonSNPs_NoNeutralVariants_CDS_2012.05.07.xlsx', sep = "/" )

  if ( file.exists( hybcappath ) ){

    require( readxl )

    data = read_excel(

      hybcappath,
      col_types = c( 
        "text",
        "numeric",
        "text",
        "numeric",
        "text",
        "numeric",
        "numeric",
        rep("text",44)
      )
    )

    dbWriteTable(
      full_con,
      "ccle_hybcap",
      as.data.frame( 
        data
      ),
      overwrite = T
    )
    
    print(paste0("Parsed file ", hybcappath))
  }

  # Cosmic CLP parsing

  cosmicclppath = paste( parser_path, 'CosmicCLP_CompleteExport.tsv', sep = "/" )
  if (file.exists( cosmicclppath) ){
    importCosmicCLP_exome(cosmicclppath, full_con)
    print(paste0("Parsed file ", cosmicclppath))
  }

  print("Parsed all available data")
}

#' Loads VCF-based data into the db
parse_vcf_data_into_db = function(){
  
}