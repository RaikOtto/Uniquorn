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
  if (file.exists( idspath) ){
    importCellLineIDs(idspath, full_con)
    print(paste0("Parsed file ", idspath))
  }

    #infopath
  infopath = paste( parser_path, 'CCLE_sample_info_file_2012-10-18.txt', sep = "/" )
  if (file.exists( infopath) ){
    importCCLE_info(infopath , full_con)
    print(paste0("Parsed file ", infopath))
  }

    # CCLE gene expression

  affypath = paste( parser_path, 'CCLE_Expression_Entrez_2012-09-29.gct', sep = "/" )
  if (file.exists( affypath) ){
    importCCLE_affy(affypath , full_con)
    print(paste0("Parsed file ", affypath))
  }

  # ccle genotype data
  importCCLE_hybcap = function (fn, con) {
    require(readxl)
    data <- read_excel(fn, col_types = c("text", "numeric", "text","numeric", "text", "numeric", "numeric", "text", "text","text", "text", "text", "text", "text", "text", "text", "text", "text", "text", "text", "text","text", "text", "text", "text", "text", "text","text", "text", "text", "text", "text", "text","text", "text", "text", "text", "text", "text","text", "text", "text", "text", "text", "text","text", "text", "text", "text", "text", "text"))
    dbWriteTable(con, "ccle_hybcap", as.data.frame(data), overwrite = TRUE)
  }
  hybcappath = paste( parser_path, 'CCLE_hybrid_capture1650_hg19_NoCommonSNPs_NoNeutralVariants_CDS_2012.05.07.xlsx', sep = "/" )

  if (file.exists( hybcappath) ){
    importCCLE_hybcap(hybcappath , full_con)
    print(paste0("Parsed file ", hybcappath))
  }

  # CCLE copy number data

  cnpath   = paste( parser_path, 'CCLE_copynumber_byGene_2013-12-03.txt', sep = "/" )
  if (file.exists( cnpath) ){
    importCCLE_cn(cnpath , full_con)
    print(paste0("Parsed file ", cnpath))
  }

  # Cosmic CLP parsing

  cosmicclppath = paste( parser_path, 'CosmicCLP_CompleteExport.tsv', sep = "/" )
  if (file.exists( cosmicclppath) ){
    importCosmicCLP_exome(cosmicclppath, full_con)
    print(paste0("Parsed file ", cosmicclppath))
  }

  # drug response

  drugpath   = paste( parser_path, 'CCLE_NP24.2009_Drug_data_2015.02.24.csv', sep = "/" )
  if (file.exists( cnpath) ){
    importCCLE_drugresponse(drugpath , full_con)
    print(paste0("Parsed file ", drugpath))
  }

  print("Parsed all available data")
}
