### parse files
#' Parses data into r list variable
#' @export
initiate_uniquorn_database = function( 
    #cosmic_genotype_file = "CosmicCLP_CompleteExport.tsv",
    cosmic_genotype_file = "CosmicCLP_MutantExport.tsv",
    cellminer_genotype_file = 'DNA__Exome_Seq_none.txt',
    ccle_genotype_file = "CCLE_hybrid_capture1650_hg19_allVariants_2012.05.07.maf",
    ucsc_db_snp_file = 'snp142Common.txt',
    ref_gen = "hg19"
  ){
  
  suppressPackageStartupMessages(library("stringr"))
  print( c( "Found CoSMIC: ", file.exists(cosmic_genotype_file) )  )
  print( c( "Found CCLE: ", file.exists(ccle_genotype_file) )  )
  print( c( "Found CellMiner: ", file.exists(cellminer_genotype_file) )  )
  print( c( "Found DbSNP: ", file.exists(ucsc_db_snp_file) )  )
  
  ### pre processing

  #db snp integration

  ref_gen_path = paste0( c( system.file("", package="Uniquorn"), "/", ref_gen ,"/" ), collapse = "")
  dir.create( ref_gen_path, showWarnings = F)

  path_to_python_dbsnp_python_parser    = paste( system.file("", package="Uniquorn"),"parse_db_snp.py", sep ="/")
  path_to_python_dbsnp_python_parser_db = paste( ref_gen_path,"parse_db_snp_python.pickle", sep ="/")
  path_to_python                        = paste( system.file("", package="Uniquorn"),"pre_compute_raw_data.py", sep ="/")
  
  if (  file.exists( ucsc_db_snp_file  ) ){
    
    print( paste0( c( "Found DbSNP file",  ucsc_db_snp_file, ", preprocessing."), collapse = " " ) )
    command_line = str_c(
      c(  
        'python', path_to_python_dbsnp_python_parser, 
        "-i", ucsc_db_snp_file,
        "-o", path_to_python_dbsnp_python_parser_db
      ), 
      collapse = " "
    )
      
    if ( ! file.exists( path_to_python_dbsnp_python_parser_db ) )
      system( command_line, ignore.stdout = F, intern = F )
    print( "Finished DbSNP pre-processing" )
  }

  # python parser
  
  command_line = str_c( 
    c(  
      'python',           path_to_python,
      "-ccle ",           ccle_genotype_file,
      "-cosmic ",         cosmic_genotype_file,
      "-cellminer",       cellminer_genotype_file,
      "-o_ref_gen_path ", ref_gen_path,
      "-i_dbsnp",         path_to_python_dbsnp_python_parser_db,
      "-filter_frequent"
    ),
    collapse = " "
  )
  
  system( command_line, ignore.stdout = F, intern = F )
  
  message("Parsing data finished")

}