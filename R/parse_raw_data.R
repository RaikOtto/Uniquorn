### parse files

#' Parses data into the sql database
#' Requires the currently (Nov 2015) github located package "chapmandu2/CancerCellLines"
#' @export
parse_data_into_younikorn_db = function( parser_path, db_path = system.file("database", package="Younikorn")  ){

  require(RSQLite)
  db_path  = paste0(db_path,"/Younikorn.db")
  print( paste0( "Parsing data and storing in db: ",db_path) )
  
  drv = dbDriver("SQLite")
  full_con = dbConnect( drv, dbname = db_path )

  # ids data
  
  idspath = system.file(
    "extdata",
    "CellLineIDNormalisationOct15.txt",
    package = "CancerCellLines"
  )
  
  if (file.exists( idspath ) ){
    
    data = read.table( idspath, header = T, sep = "\t")
    dbWriteTable( full_con, "cell_line_ids", data, overwrite = TRUE)
    print( paste0( "Parsed file ", idspath ) )
  }

  #infopath
  
  infopath = paste( parser_path, 'CCLE_sample_info_file_2012-10-18.txt', sep = "/" )
  
  if ( file.exists( infopath ) ){
    
    data = read.table( infopath, header = T, sep = "\t")
    colnames( data ) = c( 
      "CCLE_name",
      "Primary_cell_name",
      "Cell_line_aliases",
      "Gender",
      "Site_primary",
      "Histology",
      "Hist_subtype1",
      "Notes",
      "Source",
      "Expression_arrays",
      "SNP_arrays",
      "Oncomap",
      "Hybrid_capture_sequencing"
    )
    
    dbWriteTable( 
      full_con,
      "ccle_sampleinfo",
      data,
      overwrite = T
    )
    
    print( paste0( "Parsed file ", infopath ) )
  }
  
  # ccle genotype data

  hybcappath = paste( 
    parser_path,
    'CCLE_hybrid_capture1650_hg19_NoCommonSNPs_NoNeutralVariants_CDS_2012.05.07.xlsx',
    sep = "/"
  )

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
    
    print(
      paste0(
        "Parsed file ",
        hybcappath
      )
    )
  }

  # Cosmic CLP parsing

  cosmicclppath = paste( parser_path, 'CosmicCLP_CompleteExport.tsv', sep = "/" )
  
  if (file.exists( cosmicclppath) ){
    
    require(readr)
    message("Parse the Cosmic CLP exome data file")
    data = read_tsv(
      cosmicclppath,
      col_names = c(
        "gene_name",
        "accession_number",
        "hgnc_id",
        "sample_name",
        "id_sample",
        "id_tumour",
        "mutation_id",
        "mutation_cds",
        "mutation_aa",
        "mutation_description",
        "mutation_zygosity",
        "loh",
        "grch",
        "mutation_genome_position",
        "strand",
        "snp",
        "fathmm_prediction",
        "fathmm_score",
        "mutation_somatic_status"
        ),
      col_types = "cc_iccc_________cccccccccccdc________",
      skip = 1
    )
    
    message("Write the data to the database")
    
    dbWriteTable( 
      full_con,
      "cosmicclp_exome",
      as.data.frame( data ), 
      overwrite = T
    )
    
    message("Indexing the table")
    
    dbSendQuery( 
      full_con,
      sprintf(
        " CREATE INDEX `cosmicclp_exome_gene_name` ON `%s` (`gene_name` ASC); ",
        "cosmicclp_exome"
      )
    )
    
    dbSendQuery(
      full_con,
      sprintf(
        " CREATE INDEX `cosmicclp_exome_sample_name` ON `%s` (`sample_name` ASC); ",
        "cosmicclp_exome"
      )
    )
    dbSendQuery( 
      full_con,
      sprintf(
        " CREATE INDEX `cosmicclp_exome_gene_name_AND_sample_name` ON `%s` (`gene_name`,`sample_name` ASC); ",
        "cosmicclp_exome"
      )
    )
    
    #message("Finished importing Cosmic CLP exome sequencing data")
    
    print(paste0("Parsed file ", cosmicclppath))
  }

  print("Parsed all available data")
}

#' Loads VCF-based data into the db
parse_vcf_data_into_db = function(){
  
}