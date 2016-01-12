split_coords = function(vec){
    
    require("stringr", quietly = T)
    split_vec = str_split( vec, ":" )
    chrom = split_vec[1]
}

parse_cosmic_genotype_data = function( cosmic_file, sim_list ){
    
    require("stringr", quietly = T)
    
    exclude_cols_cosmic = c(rep("NULL",4),"character",rep("NULL",13),"character",rep("NULL",13))
    
    #cosmic_genotype_file = "/Users/raik_000/Dropbox/PhD/Uniquorn_project//Raw_data/CosmicCLP_MutantExport.tsv"
    cosmic_genotype_tab = read.csv2( cosmic_file, sep ="\t", colClasses = exclude_cols_cosmic)
    
    coords = as.character( sapply( cosmic_genotype_tab[,2], FUN = str_replace_all, ":|-", "_" ) )
    cls    = str_replace_all( str_to_upper(cosmic_genotype_tab[,1]), "/|(|])| ", "" )
    
    new_sim_list = data.frame( coords, cls )
    colnames(new_sim_list) = colnames(sim_list)
    sim_list = rbind( sim_list, new_sim_list )
    
    return(sim_list)
}

parse_ccle_genotype_data = function( ccle_file, sim_list ){
    
    require("stringr", quietly = T)
    
    exclude_cols_ccle = c(rep("NULL",4),rep("character", 3),rep("NULL",8),"character",rep("NULL",35))
    
    ccle_genotype_tab = read.csv2( ccle_file, sep ="\t", colClasses = exclude_cols_ccle)
    
    coords = as.character( paste0( c( ccle_genotype_tab[,c(1,2,3)] ), collapse = "_" ) )
    cls    = sapply( ccle_genotype_tab[,4], FUN = str_split, "_" )
    cls    = as.character( sapply( cls, FUN = function(vec){ return(vec[1]) }) )
    
    new_sim_list = data.frame( coords, cls )
    colnames(new_sim_list) = colnames(sim_list)
    sim_list = rbind( sim_list, new_sim_list )
    
    return(sim_list)
}

