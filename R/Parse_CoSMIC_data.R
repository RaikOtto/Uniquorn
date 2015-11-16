
#' Parse the cosmic data
#' @export
parse_cosmic_raw_data = function( parser_path ){
  
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
    
    return( data_matrix )
  }
}