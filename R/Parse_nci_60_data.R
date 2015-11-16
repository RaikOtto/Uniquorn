
#' Parse the Cellminer NCI60 dataset
#' @export
parse_cellminer_data = function( path_to_raw_data  ){
  
  raw_data = data.frame(

    "CL_ident" = character(),
    "HGNC_symbol" = character(),
    "Chr" = character(),
    "start" = character(),
    "stop" = character()
  )
  
  library( "stringr" )
  
  nci60_file = paste(
    path_to_raw_data,
    'DNA__Exome_Seq_none.txt',
    sep = "/"
  )
  
  cl_names = c("MCF7","MDA_MB_231","HS578T","BT_549","T47D","SF_268","SF_295","SF_539","SNB_19","SNB_75","U251","COLO205","HCC_2998","HCT_116","HCT_15","HT29","KM12","SW_620","CCRF_CEM","HL_60","K_562","MOLT_4","RPMI_8226","SR","LOXIMVI","MALME_3M","M14","SK_MEL_2","SK_MEL_28","SK_MEL_5","UACC_257","UACC_62","MDA_MB_435","MDA_N","A549","EKVX","HOP_62","HOP_92","NCI_H226","NCI_H23","NCI_H322M","NCI_H460","NCI_H522","IGROV1","OVCAR_3","OVCAR_4","OVCAR_5","OVCAR_8","SK_OV_3","NCI_ADR_RES","PC_3","DU_145","786_0","A498","ACHN","CAKI_1","RXF_393","SN12C","TK_10","UO_31")

  nci_data = read.table( nci60_file, header = F, skip = 11, fill = T, sep = "\t", strip.white = T)
  mut_ident = nci_data[,1]

  nci_data = nci_data[, 19: ( dim(nci_data)[2] )  ]
  colnames( nci_data ) = cl_names
  
  return( data_matrix )
}