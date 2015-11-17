
#' Parse the Cellminer NCI60 dataset
parse_cellminer_data = function( path_to_raw_data, raw_data ){
  
  nci60_file = paste(
    path_to_raw_data,
    'DNA__Exome_Seq_none.txt',
    sep = "/"
  )
  
  if ( file.exists( nci60_file )  ){

    message( paste0( "Found NCI-60 file, start parsing :", nci60_file))
    
    nci_data = read.table( 
      nci60_file, 
      header = F, 
      skip = 11, 
      fill = T, 
      sep = "\t", 
      strip.white = T, 
      stringsAsFactors = F,
      nrows = 100
    )
  
    new_nci_60_data = data.frame(
  
      "CL_ident" = character(),
      "HGNC_symbol" = character(),
      "Chr" = character(),
      "start" = character(),
      "stop" = character()
    )
    
    cl_names = c("MCF7","MDA_MB_231","HS578T","BT_549","T47D","SF_268","SF_295","SF_539","SNB_19","SNB_75","U251","COLO205","HCC_2998","HCT_116","HCT_15","HT29","KM12","SW_620","CCRF_CEM","HL_60","K_562","MOLT_4","RPMI_8226","SR","LOXIMVI","MALME_3M","M14","SK_MEL_2","SK_MEL_28","SK_MEL_5","UACC_257","UACC_62","MDA_MB_435","MDA_N","A549","EKVX","HOP_62","HOP_92","NCI_H226","NCI_H23","NCI_H322M","NCI_H460","NCI_H522","IGROV1","OVCAR_3","OVCAR_4","OVCAR_5","OVCAR_8","SK_OV_3","NCI_ADR_RES","PC_3","DU_145","786_0","A498","ACHN","CAKI_1","RXF_393","SN12C","TK_10","UO_31")
    
    extend_raw_data = function( cl_mutation_member, chromosome, mut_start, mut_stop, hgnc_symbol ){
      
      new_line = data.frame( 
        "CL_ident" = paste( cl_mutation_member, "CellMiner", sep = "_"), 
        "HGNC_symbol" = hgnc_symbol, 
        "Chr" = chromosome, 
        "start"= mut_start, 
        "stop" = mut_stop
      )
      new_nci_60_data <<- rbind( new_nci_60_data, new_line )
    }
    
    yield_information = function( vec ){
      
      chromosome  = vec[ 4 ]
      mut_start   = vec[ 5 ]
      mut_stop    = vec[ 6 ]
      hgnc_symbol = vec[ 2 ]
      
      genotype    = vec[ 19:( dim(nci_data)[ 2 ] - 1 ) ]
      member_cls  = cl_names[ genotype != "-"  ]
  
      lapply( member_cls, FUN = extend_raw_data, chromosome, mut_start, mut_stop, hgnc_symbol )
      
    }
    
    suppressMessages( apply( nci_data, MARGIN = 1, FUN = yield_information  ) )
    
    raw_data <<- rbind(
      raw_data,
      new_nci_60_data
    )
    
    message( paste0( "Parsed NCI-60 file:", nci60_file))
    
  } else {
    
    message( paste0( "Did not find NCI-60 file:", nci60_file))
  }
  
  return(raw_data)
}