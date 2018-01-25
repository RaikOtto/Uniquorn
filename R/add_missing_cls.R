#' add_missing_cls
#' 
#' @param res_table Table that contains the identification results
#' @param dif_cls Missing CLs
#' @return Results table with added missing cls
add_missing_cls = function(
    res_table,
    dif_cls
){

    for (cl in dif_cls){

        res_table = data.frame(
             "CL" = c( as.character(res_table$CL), 
               as.character(cl)),
        "CL_source" = c( as.character( res_table$CL_source),
            utils::tail(unlist(str_split(cl,"_")),1
           )),
        "Found_muts_abs" = c( as.character( res_table$Found_muts_abs ),
            as.character( "0" )),
        "Count_mutations_abs" = c( as.character( res_table$Count_mutations_abs),
                as.character( "0" )),
        "P_values"    = c( as.character( res_table$P_values),
                           as.character( "1" )),
        "Q_values"    = c( as.character( res_table$Q_values),
                           as.character( "1" )),
        "P_value_sig" = c( as.character( res_table$P_value_sig),
                           as.character( "FALSE" )),
        "Q_value_sig" = c( as.character( res_table$Q_value_sig),
                           as.character( "FALSE" ))
        )
    }
 return( res_table )
}