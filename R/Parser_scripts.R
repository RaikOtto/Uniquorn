#' Parses cosmic genotype data
#' @export
parse_cosmic_genotype_data = function( cosmic_file, sim_list ){
   
    require("stringr", quietly = TRUE, warn.conflicts = FALSE)
    
    split_coords = function(vec){
        
        require("stringr", quietly = TRUE)
        split_vec = str_split( vec, ":" )
        chrom = split_vec[1]
    }

    if (grepl("MutantExport", cosmic_file)){ # MutantExport
        
        exclude_cols_cosmic = c(rep("NULL",4),"character",rep("NULL",18),"character",rep("NULL",14))
    } else {
        warning("Warning. This is not the recommended COSMIC genotype file! The recommended file is the 'CosmicCLP_MutantExport.tsv.gz' file.")
        exclude_cols_cosmic = c(rep("NULL",4),"character",rep("NULL",13),"character",rep("NULL",13))
    }
    
    cosmic_genotype_tab = read.csv2( cosmic_file, sep ="\t", colClasses = exclude_cols_cosmic)
    
    coords = as.character( sapply( cosmic_genotype_tab[,2], FUN = str_replace_all, ":|-", "_" ) )
    cls    = str_replace_all( str_to_upper(cosmic_genotype_tab[,1]), "/|(|])| ", "" )
    cls    = sapply( cls, FUN = function( cl_name ){ return( paste(cl_name, "COSMIC", sep = "_") ) } )
    
    new_sim_list = data.frame( coords, cls )
    colnames(new_sim_list) = colnames(sim_list)
    sim_list = rbind( sim_list, new_sim_list )
    
    return(sim_list)
}

#' Parses ccle genotype data
#' @export
parse_ccle_genotype_data = function( ccle_file, sim_list ){
    
    require("stringr", quietly = T)
    
    exclude_cols_ccle = c(rep("NULL",4),rep("character", 3),rep("NULL",8),"character",rep("NULL",35))
    
    ccle_genotype_tab = read.csv2( ccle_file, sep ="\t", colClasses = exclude_cols_ccle)
    
    coords = as.character( apply(
        ccle_genotype_tab,
        FUN = function( vec ){ return( paste0( c(
            as.character( vec[1] ),
            as.character( vec[2] ),
            as.character( vec[3] )
            ), collapse = "_" ) ) },
        MARGIN = 1
    ) )
    
    cls    = sapply( ccle_genotype_tab[,4], FUN = str_split, "_" )
    cls    = as.character( sapply( cls, FUN = function(vec){ return(vec[1]) }) )
    cls    = as.character( sapply( cls, FUN = function( cl_name ){ return( paste(cl_name, "CCLE", sep = "_") ) } ) )
    
    new_sim_list = data.frame( coords, cls )
    colnames(new_sim_list) = colnames(sim_list)
    sim_list = rbind( sim_list, new_sim_list )
    
    return(sim_list)
}

